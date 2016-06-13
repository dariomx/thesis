import numpy as np
from test_util import conv_mat

# http://computationalmathematics.org/topics/files/IC.html
def ichol(M):
    A = M.toarray()
    n = A.shape[0]
    for k in xrange(n):
        A[k,k] = np.sqrt(A[k,k])
        for i in xrange(k+1,n):
            if A[i,k] != 0:
                A[i,k] = A[i,k] / A[k,k]
        for j in xrange(k+1,n):
            for i in xrange(j,n):
                if A[i,j] != 0:
                    A[i,j] = A[i,j] - A[i,k]*A[j,k]
    return conv_mat(A, M.getfmt())

