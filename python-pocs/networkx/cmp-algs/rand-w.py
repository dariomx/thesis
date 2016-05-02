#!/usr/bin/env python

from sys import argv
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import scipy.stats
from test_util import eprint

def rand_w(dist_name, dist_params, n, pz, fn):
    print("generating random weight matrix with dist %s" % dist_name)
    dist = getattr(scipy.stats, dist_name)
    dps = map(float, dist_params.split(","))
    tril = np.tril_indices(n, -1)
    m = len(tril[0])
    R = dist.rvs(*dps[:-2], loc=dps[-2], scale=dps[-1], size=m)
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
    if len(argv) < 5:
        args = "<pz> <dist_name> <dist_params> <n1> [<n2> <n3>...]"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    pz = float(argv[1])
    dist_name = argv[2]
    dist_params = argv[3]
    ns  = map(int, argv[4:])
    fns = map(lambda n: str(n) + ".mtx", ns)
    for n,fn in zip(ns, fns):
        rand_w(dist_name, dist_params, n, pz, fn)
