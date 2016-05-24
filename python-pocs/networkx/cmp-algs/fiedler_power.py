from scipy.linalg import norm
from scipy.sparse.linalg import spsolve as solve
from scipy.sparse import eye
from numpy import dot
from test_util import get_rand_dist, relres, eprint
from numpy import ones, loadtxt

# computes fiedler eigenpair using rqi (provides initial estimations)
def fiedler_rqi(L, tol=1e-7):
    dist = get_rand_dist("uniform","0,1")
    ac = 0.5
    #fv = dist(L.shape[0])
    fv = loadtxt("domain/fiedler-vec.mat") + dist(L.shape[0])
    return rqi(L, ac, fv, tol)

def power_method(A, v, tol):
    rr = 1
    while rr > tol:
        v = v / norm(v)
        v = A * v
        lam = dot(v, A*v)
        rr = relres(A, v, lam)
        eprint("lam=%.3E, rr=%.3E" % (lam, rr))
    return lam, v

def fiedler_pm(L, tol=1e-7):
    dist = get_rand_dist("uniform","0,1")
    fv = dist(L.shape[0])
    return power(method(L, fv))
