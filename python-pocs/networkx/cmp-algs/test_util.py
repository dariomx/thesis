from __future__ import print_function
from sys import stderr
from datetime import datetime
from itertools import izip
from numpy import loadtxt, sign, ndarray, concatenate, inf
from numpy import diagflat, diagonal, tril_indices
from numpy.linalg import cond
from scipy.linalg import norm as dnorm, eigh, orth
from scipy.sparse.linalg import norm as snorm
from scipy.sparse import issparse, lil_matrix, csr_matrix
import scipy.stats
import numpy as np
from scipy.io import mmread, mminfo, savemat, mmwrite
import networkx as nx
import matplotlib.pyplot as plt

mat_norm_ord = 'fro'
vec_norm_ord = 2
    
def eprint(*args, **kwargs):
    print(*args, file=stderr, **kwargs)

def conv_mat(M, fmt):
    if fmt == "dense":
        if type(M) == ndarray:
            Mc = M
        else:
            Mc = M.toarray()
    else:
        import scipy
        method = getattr(scipy.sparse, fmt + "_matrix")
        Mc = method(M)
    return Mc

def load_mat(fn, fmt):
    eprint("loading matrix ... ")
    start = datetime.now()
    if fn.endswith(".mat"):
        M = loadtxt(fn)
    elif fn.endswith(".mtx") or fn.endswith(".mtz.gz"):
        eprint("mm header:  " + str(mminfo(fn)))
        M = mmread(fn)
    end = datetime.now()
    load_time = end - start
    eprint("loaded in %s" % load_time)
    start = datetime.now()
    Mc = conv_mat(M, fmt)
    end = datetime.now()
    conv_time = end - start
    eprint("converted to format %s in %s" % (fmt, conv_time))
    return Mc    

def load_fiedler_vec(ac_fn, fv_fn):
    ac = None if ac_fn is None else loadtxt(ac_fn)[1]
    fv = None if fv_fn is None else loadtxt(fv_fn)
    return ac, fv

def is_nearly_zero(x):
    eps = 2.220446049250313E-16
    zero_tol = 100000.0 * eps
    return True if abs(x) < zero_tol else False

# loads a graph from its laplacian
def load_graph(fn, is_lap, zeros):
    M = load_mat(fn, fmt="dense")
    G = nx.Graph()
    G.add_nodes_from(xrange(1, M.shape[0]+1))
    Mc = M.tocoo()
    for i, j, lw in izip(Mc.row, Mc.col, Mc.data):
        if is_lap:
            w = 0 if i==j else -lw
        else:
            w = lw
        zw = 0 if zeros and is_nearly_zero(w) else w
        G.add_edge(i, j, weight=zw)
    return G

# invert y if the signs are opposite as x
def invsign(y, x):
    return -y if (sign(x) == -sign(y)).all() else y

def relres(L, ac, fv, Lnorm=None):
    norm = snorm if issparse(L) else dnorm
    Lnorm = norm(L, mat_norm_ord) if Lnorm is None else Lnorm
    return dnorm(L.dot(fv) - ac*fv, vec_norm_ord) / Lnorm
        
def relerr_vec(x, y):
    return dnorm(x - invsign(y, x)) / dnorm(x)

def relerr_mat(A, B):
    assert issparse(A) == issparse(B), "incompatible matrices"
    norm = snorm if issparse(A) else dnorm
    return norm(A - B, mat_norm_ord) / norm(A, mat_norm_ord)

def cmp_ac_fv(L, cac, cfv, ac, fv):
    eprint("alg conn: %.16f" % cac)
    eprint("relres c: %.16f" % relres(L, cac, cfv))
    if ac is not None and fv is not None:
        eprint("relres i: %.16f" % relres(L, ac, fv))
        args = (relerr(ac, cac), relerr(fv, cfv))
        eprint("relerr:  ac = %.16f,  fv = %.16f" % args)

# calculates the (unormalized) laplacian from the weights matrix
def lap(W, fmt):
    d = W.sum(axis=1)
    d = d if len(d.shape) == 1 else concatenate(d.A)
    D = lil_matrix(W.shape)
    D.setdiag(d)
    L = csr_matrix(D) - csr_matrix(W)
    return conv_mat(L, fmt)

# upper bound of ac discovered by Fiedler, assumes similary (weight) graph 
# had 1's on the diagonal (so they got lost while doing L = D - W, and
# we recover here: diag(D) = diag(L) + diag(W) = diag(L) + diag(I))
def get_ac_upbound(L):
    n = float(L.shape[0])
    dmin = min(map(lambda x: x+1, L.diagonal()))
    return n/(n-1), dmin, n/(n-1) * dmin

# returns the original weights from the laplacian (just the lower half
# given symmetry + plus the diagonal).
def get_weights(L):
    n = L.shape[0]
    ix, iy = tril_indices(n, 0)
    return ix, iy, -L[ix, iy]

def get_lap(fn, fmt):
    W = load_mat(fn, fmt)
    if is_lap(fn):
        L = W
    else:
        L = lap(W, fmt)
    return L

def cond_err(A):
    ls, V = eigh(conv_mat(A, fmt="dense"), overwrite_a=False)
    D = diagflat(ls)
    return cond(V, p=mat_norm_ord), relerr_mat(A, V.T*D*V)

def get_eigvals(L):
    n = L.shape[0]
    Ld = conv_mat(L, fmt="dense")
    ls, _ = eigh(L, eigvals=(0,n-1), overwrite_a=True)
    return ls

# returns a distribution object and its parsed params
def get_dist(dist_name, dist_params):
    print("random generator will use dist %s" % dist_name)
    dist = getattr(scipy.stats, dist_name)
    dps = map(float, dist_params.split(","))
    return dist, dps

# returns a random generator for given distribution/params
def get_rand_dist(dist_name, dist_params):
    dist, dps = get_dist(dist_name, dist_params)
    return lambda s: dist.rvs(*dps[:-2], loc=dps[-2], scale=dps[-1], size=s)

# returns the probability density function of a distribution/params
def get_pdf_dist(dist_name, dist_params):
    dist, dps = get_dist(dist_name, dist_params)
    return lambda xs: dist.pdf(xs, *dps[:-2], loc=dps[-2], scale=dps[-1])

# returns random matrix with given eigenvalues
# assumes that rand_dist produces non-singular matrices
def rand_mat_eigv(rand_dist, n, eigv):
    D = np.diag(eigv)
    Q = orth(rand_dist((n,n)))
    M = np.dot(Q, np.dot(D, Q.T))
    return csr_matrix(M)

# tells if the matrix represents a laplacian base on file name convention
def is_lap(fn):
    return fn.endswith("-lap.mtx")

def load_weights(fn):
    L = get_lap(fn, fmt="csr")
    eprint("extracting weights from %s ... " % fn)        
    ix, iy, ws = get_weights(L)
    nnz = np.squeeze(np.asarray(ws[ws > 0].T))
    ixy = izip(ix, iy)    
    ix_nnz, iy_nnz = izip(*[ii for (ii,w) in izip(ixy, nnz) if w > 0])
    args = (ws.shape[1], nnz.shape[0])
    eprint("extracted %d weights (%d nnz)" % args)
    return ix_nnz, iy_nnz, nnz.shape[0], nnz

def get_plot_axis(is_3d=False):
    fig = plt.figure()
    if is_3d:
        return fig.add_subplot(111, projection="3d")
    else:
        return fig.add_subplot(111)

def save_plot(ax, title, img_file, bgcolor="black", fgcolor="white"):
    ax.set_axis_bgcolor(bgcolor)
    ax.spines['bottom'].set_color(fgcolor)
    ax.spines['top'].set_color(fgcolor)
    ax.spines['right'].set_color(fgcolor)
    ax.spines['left'].set_color(fgcolor)
    ax.tick_params(axis='x', colors=fgcolor)
    ax.tick_params(axis='y', colors=fgcolor)
    ax.yaxis.label.set_color(fgcolor)
    ax.xaxis.label.set_color(fgcolor)
    ax.title.set_color(fgcolor)
    ax.set_title(title)
    plt.savefig(img_file, facecolor=bgcolor, edgecolor='none')

parse_bool = lambda s: s == "true" or s == "True"
