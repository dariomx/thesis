from __future__ import print_function
from sys import stderr
from datetime import datetime
from itertools import izip
from numpy import loadtxt, sign, ndarray, concatenate, inf
from numpy import diagflat
from numpy.linalg import cond
from scipy.linalg import norm as dnorm, eigh
from scipy.sparse.linalg import norm as snorm
from scipy.sparse import issparse, lil_matrix, csr_matrix
from scipy.io import mmread, mminfo, savemat
import networkx as nx

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

def relres(L, ac, fv):
    rr = None
    if issparse(L):
        rr = dnorm(L*fv - ac*fv, vec_norm_ord) / snorm(L, mat_norm_ord)
    else:
        rr = dnorm(L.dot(fv) - ac*fv, vec_norm_ord) / dnorm(L, mat_norm_ord)
    return rr
        
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
    
def lap(W, fmt):
    d = W.sum(axis=1)
    d = d if len(d.shape) == 1 else concatenate(d.A)
    D = lil_matrix(W.shape)
    D.setdiag(d)
    L = csr_matrix(D) - csr_matrix(W)
    return conv_mat(L, fmt)

def get_lap(fn, is_lap, fmt):
    if is_lap:
        L = load_mat(fn, fmt)
    else:
        W = load_mat(fn, fmt)
        L = lap(W, fmt)
    return L

def cond_err(A):
    ls, V = eigh(conv_mat(A, fmt="dense"), overwrite_a=False)
    D = diagflat(ls)
    return cond(V, p=mat_norm_ord), relerr_mat(A, V.T*D*V)

parse_bool = lambda s: s == "true" or s == "True"
