from datetime import datetime
from itertools import izip
from numpy import loadtxt, sign, ndarray
from scipy.linalg import norm as dnorm
from scipy.sparse.linalg import norm as snorm
from scipy.sparse import issparse
from scipy.io import mmread, mminfo
import networkx as nx

def load(L_file, dense, ac_file=None, fv_file=None):
    if L_file.endswith(".mat"):
        L  = loadtxt(L_file)
    elif L_file.endswith(".mtx") or L_file.endswith(".mtz.gz"):
        print "loading mm format:  " + str(mminfo(L_file))
        start = datetime.now()
        L = mmread(L_file)
        L = L.toarray() if dense else L.tocsr() 
        end = datetime.now()
        print "loaded in " + str(end-start)
    ac = None if ac_file is None else loadtxt(ac_file)[1]
    fv = None if fv_file is None else loadtxt(fv_file)
    return L, ac, fv

def is_nearly_zero(x):
    eps = 2.220446049250313E-16
    zero_tol = 100000.0 * eps
    return True if abs(x) < zero_tol else False

# loads a graph from its laplacian
def load_graph(L_file, zeros=False):
    L, _, _ = load(L_file, dense=False)
    G = nx.Graph()
    G.add_nodes_from(xrange(1, L.shape[0]+1))
    Lc = L.tocoo()
    for i, j, w in izip(Lc.row, Lc.col, Lc.data):
        zw = 0 if zeros and is_nearly_zero(w) else w
        G.add_edge(i, j, weight=zw)
    return G

# invert y if the signs are opposite as x
def invsign(y, x):
    return -y if (sign(x) == -sign(y)).all() else y

def relres(L, ac, fv):
    norm = snorm if issparse(L) else dnorm
    return dnorm(L*fv - ac*fv) / norm(L)

def relerr(x, y):
    return dnorm(x - invsign(y, x)) / dnorm(x)

def cmp(L, ac, fv, cac, cfv):
    print "alg conn: %.16f" % cac
    print "relres c: %.16f" % relres(L, cac, cfv)
    if ac is not None and fv is not None:
        print "relres i: %.16f" % relres(L, ac, fv)
        args = (relerr(ac, cac), relerr(fv, cfv))
        print "relerr:  ac = %.16f,  fv = %.16f" % args
    
def test(argv, calc, dense=False):
    if len(argv) == 4:
        L, ac, fv = load(argv[1], dense, argv[2], argv[3])
    else:
        L, ac, fv = load(argv[1], dense)  
    start = datetime.now()
    cac, cfv = calc(L)
    end = datetime.now()
    print "calc took %s" % (end - start)
    cmp(L, ac, fv, cac, cfv)
