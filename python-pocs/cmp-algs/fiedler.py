from re import compile
from functools import partial
from numpy import (array, asmatrix, asarray, dot, matrix, ndarray, ones,
                   reshape, sqrt, zeros, random)
from numpy.linalg import norm, qr
from numpy.random import normal
from scipy.linalg import eigh, inv
from scipy.sparse import csc_matrix, spdiags
from scipy.sparse.linalg import eigsh, lobpcg
from scipy.linalg.blas import dasum, ddot, daxpy

_tracemin_method = compile('^tracemin(?:_(.*))?$')

class _PCGSolver(object):
    """Preconditioned conjugate gradient method.
    """

    def __init__(self, A, M):
        self._A = A
        self._M = M or (lambda x: x.copy())

    def solve(self, B, tol):
        B = asarray(B)
        X = ndarray(B.shape, order='F')
        for j in range(B.shape[1]):
            X[:, j] = self._solve(B[:, j], tol)
        return X

    def _solve(self, b, tol):
        A = self._A
        M = self._M
        tol *= dasum(b)
        # Initialize.
        x = zeros(b.shape)
        r = b.copy()
        z = M(r)
        rz = ddot(r, z)
        p = z.copy()
        # Iterate.
        while True:
            Ap = A(p)
            alpha = rz / ddot(p, Ap)
            x = daxpy(p, x, a=alpha)
            r = daxpy(Ap, r, a=-alpha)
            if dasum(r) < tol:
                return x
            z = M(r)
            beta = ddot(r, z)
            beta, rz = beta / rz, beta
            p = daxpy(p, z, a=beta)


class _CholeskySolver(object):
    """Cholesky factorization.
    """

    def __init__(self, A):
        if not self._cholesky:
            raise ValueError('Cholesky solver unavailable.')
        self._chol = self._cholesky(A)

    def solve(self, B):
        return self._chol(B)

    try:
        from scikits.sparse.cholmod import cholesky
        _cholesky = cholesky
    except ImportError:
        _cholesky = None


class _LUSolver(object):
    """LU factorization.
    """

    def __init__(self, A):
        if not self._splu:
            raise ValueError('LU solver unavailable.')
        self._LU = self._splu(A)

    def solve(self, B):
        B = asarray(B)
        X = ndarray(B.shape, order='F')
        for j in range(B.shape[1]):
            X[:, j] = self._LU.solve(B[:, j])
        return X

    try:
        from scipy.sparse.linalg import splu
        _splu = partial(splu, permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0.,
                        options={'Equil': True, 'SymmetricMode': True})
    except ImportError:
        _splu = None

def _tracemin_fiedler(L, X, normalized, tol, method):
    """Compute the Fiedler vector of L using the TraceMIN-Fiedler algorithm.
    """
    n = X.shape[0]

    if normalized:
        # Form the normalized Laplacian matrix and determine the eigenvector of
        # its nullspace.
        e = sqrt(L.diagonal())
        D = spdiags(1. / e, [0], n, n, format='csr')
        L = D * L * D
        e *= 1. / norm(e, 2)

    if not normalized:
        def project(X):
            """Make X orthogonal to the nullspace of L.
            """
            X = asarray(X)
            for j in range(X.shape[1]):
                X[:, j] -= X[:, j].sum() / n
    else:
        def project(X):
            """Make X orthogonal to the nullspace of L.
            """
            X = asarray(X)
            for j in range(X.shape[1]):
                X[:, j] -= dot(X[:, j], e) * e


    if method is None:
        method = 'pcg'
    if method == 'pcg':
        # See comments below for the semantics of P and D.
        def P(x):
            x -= asarray(x * X * X.T)[0, :]
            if not normalized:
                x -= x.sum() / n
            else:
                x = daxpy(e, x, a=-ddot(x, e))
            return x
        solver = _PCGSolver(lambda x: P(L * P(x)), lambda x: D * x)
    elif method == 'chol' or method == 'lu':
        # Convert A to CSC to suppress SparseEfficiencyWarning.
        A = csc_matrix(L, dtype=float, copy=True)
        # Force A to be nonsingular. Since A is the Laplacian matrix of a
        # connected graph, its rank deficiency is one, and thus one diagonal
        # element needs to modified. Changing to infinity forces a zero in the
        # corresponding element in the solution.
        i = (A.indptr[1:] - A.indptr[:-1]).argmax()
        A[i, i] = float('inf')
        solver = (_CholeskySolver if method == 'chol' else _LUSolver)(A)
    else:
        raise ValueError('unknown linear system solver.')

    # Initialize.
    Lnorm = abs(L).sum(axis=1).flatten().max()
    project(X)
    W = asmatrix(ndarray(X.shape, order='F'))

    while True:
        # Orthonormalize X.
        X = qr(X)[0]
        # Compute interation matrix H.
        W[:, :] = L * X
        H = X.T * W
        sigma, Y = eigh(H, overwrite_a=True)
        # Compute the Ritz vectors.
        X *= Y
        # Test for convergence exploiting the fact that L * X == W * Y.
        res = dasum(W * asmatrix(Y)[:, 0] - sigma[0] * X[:, 0]) / Lnorm
        if res < tol:
            break
        # Depending on the linear solver to be used, two mathematically
        # equivalent formulations are used.
        if method == 'pcg':
            # Compute X = X - (P * L * P) \ (P * L * X) where
            # P = I - [e X] * [e X]' is a projection onto the orthogonal
            # complement of [e X].
            W *= Y  # L * X == W * Y
            W -= (W.T * X * X.T).T
            project(W)
            # Compute the diagonal of P * L * P as a Jacobi preconditioner.
            D = L.diagonal()
            D += 2. * (asarray(X) * asarray(W)).sum(axis=1)
            D += (asarray(X) * asarray(X * (W.T * X))).sum(axis=1)
            D[D < tol * Lnorm] = 1.
            D = 1. / D
            # Since TraceMIN is globally convergent, the relative residual can
            # be loose.
            X -= solver.solve(W, 0.1)
        else:
            # Compute X = L \ X / (X' * (L \ X)). L \ X can have an arbitrary
            # projection on the nullspace of L, which will be eliminated.
            W[:, :] = solver.solve(X)
            project(W)
            X = (inv(W.T * X) * W.T).T  # Preserves Fortran storage order.

    return sigma, asarray(X)


def _get_fiedler_func(method):
    """Return a function that solves the Fiedler eigenvalue problem.
    """
    match = _tracemin_method.match(method)
    if match:
        method = match.group(1)
        def find_fiedler(L, x, normalized, tol):
            q = 2 if method == 'pcg' else min(4, L.shape[0] - 1)
            X = asmatrix(normal(size=(q, L.shape[0]))).T
            sigma, X = _tracemin_fiedler(L, X, normalized, tol, method)
            return sigma[0], X[:, 0]
    elif method == 'lanczos' or method == 'lobpcg':
        def find_fiedler(L, x, normalized, tol):
            L = csc_matrix(L, dtype=float)
            n = L.shape[0]
            if normalized:
                D = spdiags(1. / sqrt(L.diagonal()), [0], n, n, format='csc')
                L = D * L * D
            if method == 'lanczos' or n < 10:
                # Avoid LOBPCG when n < 10 due to
                # https://github.com/scipy/scipy/issues/3592
                # https://github.com/scipy/scipy/pull/3594
                sigma, X = eigsh(L, 2, which='SM', tol=tol,
                                 return_eigenvectors=True)
                return sigma[1], X[:, 1]
            else:
                X = asarray(asmatrix(x).T)
                M = spdiags(1. / L.diagonal(), [0], n, n)
                Y = ones(n)
                if normalized:
                    Y /= D.diagonal()
                sigma, X = lobpcg(L, X, M=M, Y=asmatrix(Y).T, tol=tol,
                                  maxiter=n, largest=False)
                return sigma[0], X[:, 0]
    else:
        raise ValueError("unknown method '%s'." % method)

    return find_fiedler

def fiedler_vector(L, normalized=False, tol=1e-8, method='tracemin'):
    find_fiedler = _get_fiedler_func(method)
    x = None if method != 'lobpcg' else random.rand(L.shape[0])
    return find_fiedler(L, x, normalized, tol)
