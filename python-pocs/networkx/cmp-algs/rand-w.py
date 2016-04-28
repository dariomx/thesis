#!/usr/bin/env python

from sys import argv
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmwrite

def rand_w(n, pz, fn):
    print "generating random weight matrix ... "
    tril = np.tril_indices(n, -1)
    m = len(tril[0])
    R = np.random.random(m)
    R[0:int(m*pz)] = 0
    for i in xrange(3):
        np.random.shuffle(R)
    W = np.zeros((n,n))
    W[tril] = R
    W = W + W.T
    np.fill_diagonal(W, 1)
    print "saving matrix to " + fn
    mmwrite(fn, csr_matrix(W))

# main
if __name__ == '__main__':
    pz = float(argv[1])
    ns  = map(int, argv[2:])
    fns = map(lambda n: str(n) + ".mtx", ns)
    for n,fn in zip(ns, fns):
        rand_w(n, pz, fn)
