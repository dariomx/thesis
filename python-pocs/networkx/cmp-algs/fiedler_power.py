from scipy.linalg import norm
from scipy.sparse.linalg import factorized, spsolve
from scipy.linalg import lu_factor, lu_solve, solve
from scipy.sparse import eye, issparse, csc_matrix
from numpy import dot, ones
from test_util import get_rand_dist, relres, eprint, mat_norm, get_ac_upbound, take_time
from numpy import loadtxt
from scipy.linalg import cho_factor, cho_solve
from sksparse.cholmod import cholesky
from scipy.sparse.linalg.interface import LinearOperator
import numpy as np
 
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
    return invpow(L, fv, get_lu_solver1(L), tol)

def fiedler_ship(L, tol=1e-7):
    n = L.shape[0]
    a = 1
    L1 = L - a*eye(n)
    fv = rand_vec(n)    
    return invpow(L1, fv, get_lu_solver1(L1), tol)

def fiedler_suip(L, tol=1e-7):
    L1 = get_spec_upd(L)
    fv = rand_vec(L.shape[0])
    return invpow(L1, fv, get_lu_solver1(L1), tol)

def fiedler_rsuip(L, tol=1e-3):
    L1 = get_spec_upd(L)
    solver = get_lu_solver1(L1)
    fv = rand_vec(L.shape[0])
    for i in xrange(4):
        ac, fv = invpow(L1, fv, solver, tol)
        tol *= 1e-1        
    return ac, fv

def fiedler_suipc(L, tol=1e-7):
    a, b, v1 = get_spec_upd(L, sep=True)
    fv = rand_vec(L.shape[0])
    return invpow(L, fv, get_chol_solver(L, a, b, v1), tol)

def rqi(A, v, solve, tol):
    i = 0
    rr = 1
    normA = mat_norm(A)
    lam = 0
    n = A.shape[0]
    while rr > tol:
        i = i + 1
        v = v / norm(v)        
        sAiv = solve(A - lam*eye(n), v)
        mu = dot(v, sAiv)
        lam = lam + 1 / mu
        rr = relres(A, lam, v, normA)
        eprint("iter=%d, lam=%.14f, rr=%.14f" % (i, lam, rr))
        v = sAiv
    return lam, v

def fiedler_surqi(L, tol=1e-7):
    L1 = get_spec_upd(L)
    fv = rand_vec(L.shape[0])
    return rqi(L1, fv, get_lu_solver2(L1), tol)

def rand_vec(n):
    dist = get_rand_dist("uniform","0,1")
    return dist(n)

def get_lu_solver1(A):
    if issparse(A):
        return factorized(A)
    else:
        lu_piv = lu_factor(A)
        return lambda b: lu_solve(lu_piv, b)

def get_lu_solver2(A):
    if issparse(A):
        return spsolve
    else:
        return solve

def get_chol_solver(A, a, b, v):
    if issparse(A):
        return get_chol_suops(A, a, b, v).solver
    else:
        raise ValueError("dense matrix not supported for now")
    
class MatSolverOp(LinearOperator):
    def __init__(self, A, solver):
        self.shape = A.shape
        self.dtype = A.dtype
        self.isreal = not np.issubdtype(self.dtype, np.complexfloating)
        self.solver = solver
        self.solve_iter = 0
        self.solve_time = 0

    def _matvec(self, x):
        solve = lambda: self.solver(x)
        result, time = take_time(solve)
        self.solve_iter += 1
        self.solve_time += time
        return result

def get_lu_op(A):
    eprint("A is sparse? %s" % (issparse(A)))
    fact = lambda: get_lu_solver1(A)
    solver, time = take_time(fact)
    eprint("lu factorization took %10.8f" % time)
    return MatSolverOp(A, solver)

def get_chol_opd(A):
    eprint("A is sparse? %s" % (issparse(A)))
    fact = lambda: cho_factor(A)
    cholf, time = take_time(fact)
    eprint("cholesky factorization took %10.8f" % time)
    return MatSolverOp(A, lambda b: cho_solve(cholf, b))

def get_chol_ops(A, a):
    eprint("A is sparse? %s" % (issparse(A)))
    fact = lambda: cholesky(A, beta=a)
    solver, time = take_time(fact)
    eprint("cholesky factorization took %10.8f" % time)
    return MatSolverOp(A, solver)

def get_chol_suops(A, a, b, v):
    eprint("A is sparse? %s" % (issparse(A)))
    fact = lambda: cholesky(A, beta=a)
    solver, time = take_time(fact)
    eprint("cholesky factorization took %10.8f" % time)
    def upd():
        V = csc_matrix(b * dot(v, v.T))
        solver.update_inplace(V)
        return None
    _, time = take_time(upd)
    eprint("spectral update took %10.8f" % time)
    return MatSolverOp(A, solver)

def get_spec_upd(L, c=1, sep=False):
    n = L.shape[0]
    a = 1e-2
    b = (get_ac_upbound(L, a) + c)/n
    v1 = ones(n)[None].T    
    if sep:
        return a, b, v1
    else:
        La = L + a*eye(n)
        V1 = dot(v1, v1.T)        
        return La + b*V1, a

# above was just for curiosity, it performs badly
# https://www.quantstart.com/articles/Cholesky-Decomposition-in-Python-and-NumPy
def spa_cho_factor(A):
    n = A.shape[0]
    L = np.zeros((n,n))
    for i in xrange(n):
        eprint("row %d" % i)
        for k in xrange(i+1):
            eprint("cell %d,%d" % (i,k))            
            tmp_sum = np.sum(L[i,j] * L[k,j] for j in xrange(k))
            if (i == k):
                L[i, k] = np.sqrt(A[i,i] - tmp_sum)
            else:
                L[i, k] = (1.0 / L[k,k] * (A[i,k] - tmp_sum))
    return csc_matrix(L)
