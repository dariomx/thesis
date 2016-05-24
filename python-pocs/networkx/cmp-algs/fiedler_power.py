from scipy.linalg import norm
from scipy.sparse.linalg import factorized
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse import eye, issparse
from numpy import dot
from test_util import get_rand_dist, relres, eprint, mat_norm
from numpy import ones, loadtxt

def power_method(A, v, tol):
    i = 0
    rr = 1
    normA = mat_norm(A)
    while rr > tol:
        i = i + 1
        v = v / norm(v)
        Av = A * v
        lam = dot(v, Av)
        rr = relres(A, lam, v, normA)
        eprint("iter=%d, lam=%.14f, rr=%.14f" % (i, lam, rr))
        v = Av        
    return lam, v

def fiedler_pm(L, tol=1e-7):
    fv = rand_vec(L.shape[0])
    return power_method(L, fv, tol)

def invpow(A, v, solve, tol):
    i = 0
    rr = 1
    normA = mat_norm(A)
    while rr > tol:
        i = i + 1
        v = v / norm(v)
        Aiv = solve(v)
        mu = dot(v, Aiv)
        lam = 1 / mu
        rr = relres(A, lam, v, normA)
        eprint("iter=%d, lam=%.14f, rr=%.14f" % (i, lam, rr))
        v = Aiv        
    return lam, v

def fiedler_invpow(L, tol=1e-7):
    fv = rand_vec(L.shape[0])
    return invpow(L, fv, get_solver(L), tol)

def fiedler_stip(L, tol=1e-7):
    n = L.shape[0]
    v1 = ones(n)[None].T
    M1 = dot(v1, v1.T)
    s1 = 20
    L1 = L - s1*M1
    fv = rand_vec(n)
    return invpow(L1, fv, get_solver(L1), tol)

def fiedler_ship(L, tol=1e-7):
    n = L.shape[0]
    sig = 10
    L1 = L - sig*eye(n)
    fv = rand_vec(n)    
    return invpow(L1, fv, get_solver(L1), tol)

def get_solver(A):
    if issparse(A):
        return factorized(A)
    else:
        lu_piv = lu_factor(A)
        return lambda b: lu_solve(lu_piv, b)

def rand_vec(n):
    dist = get_rand_dist("uniform","0,1")
    return dist(n)
    
