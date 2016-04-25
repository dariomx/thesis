#!/usr/bin/env python

from sys import argv
from random import random
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmwrite

def rand_w(n, pz, W_file):
    print "generating random weight matrix ... "
    W = np.random.random((n, n))
    tril = np.tril_indices(n, -1)
    rand_z = lambda x: 0 if random() <= pz else x
    W[tril] = map(rand_z, W[tril])
    W[np.triu_indices(n, 1)] = 0
    W = W + W.T
    np.fill_diagonal(W, 1)
    print "saving matrix to " + W_file
    mmwrite(W_file, csr_matrix(W))

# main
if __name__ == '__main__':
    n  = int(argv[1])
    pz = float(argv[2])
    W_file = argv[3]
    rand_w(n, pz, W_file)
