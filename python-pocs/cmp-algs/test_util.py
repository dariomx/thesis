from datetime import datetime
from itertools import izip
from numpy import loadtxt, sign, ndarray, concatenate, inf
from numpy import diag
from scipy.linalg import norm as dnorm
from scipy.sparse.linalg import norm as snorm
from scipy.sparse import issparse, coo_matrix, csr_matrix, lil_matrix
from scipy.io import mmread, mminfo, savemat
import networkx as nx

def load_mat(fn, dense):
    if fn.endswith(".mat"):
        M = loadtxt(fn)
        M = M if dense else csr_matrix(M)
    elif fn.endswith(".mtx") or fn.endswith(".mtz.gz"):
        print "loading mm format:  " + str(mminfo(fn))
        start = datetime.now()
        M = mmread(fn)
        M = M.toarray() if dense else M.tocsr()
        end = datetime.now()
        print "loaded in " + str(end-start)
    return M

def load_lap(L_file, dense, ac_file=None, fv_file=None):
    L  = load_mat(L_file, dense)
    ac = None if ac_file is None else loadtxt(ac_file)[1]
    fv = None if fv_file is None else loadtxt(fv_file)
    return L, ac, fv

def is_nearly_zero(x):
    eps = 2.220446049250313E-16
    zero_tol = 100000.0 * eps
    return True if abs(x) < zero_tol else False

# loads a graph from its laplacian
def load_graph(M_file, is_lap, zeros):
    M = load_mat(M_file, dense=False)
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
        rr = dnorm(L*fv - ac*fv, inf) / snorm(L, inf)
    else:
        rr = dnorm(L.dot(fv) - ac*fv, inf) / dnorm(L, inf)        
    return rr
        
def relerr(x, y):
    return dnorm(x - invsign(y, x)) / dnorm(x)

def cmp_ac_fv(L, cac, cfv, ac, fv):
    print "alg conn: %.16f" % cac
    print "relres c: %.16f" % relres(L, cac, cfv)
    if ac is not None and fv is not None:
        print "relres i: %.16f" % relres(L, ac, fv)
        args = (relerr(ac, cac), relerr(fv, cfv))
        print "relerr:  ac = %.16f,  fv = %.16f" % args
    
def lap(W, dense):
    d = W.sum(axis=1)
    d = d if len(d.shape) == 1 else concatenate(d.A)
    D = lil_matrix(W.shape)
    D.setdiag(d)
    L = csr_matrix(D) - csr_matrix(W)
    return L.toarray() if dense else L

parse_bool = lambda s: s == "true" or s == "True"
