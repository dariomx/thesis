#!/usr/bin/env python

from sys import argv
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import scipy.stats
from test_util import eprint

# ensures an upperbound of the algebraic connectivity
# (may break the contract with self-similarity)
def ensure_upbound_ac(W, upbound_ac):
    n = W.shape[0]
    i = np.random.randint(0, n)
    j = np.random.randint(0, n)
    W[i,:] = 0
    W[:,i] = 0
    W[i,j] = upbound_ac
    W[j,i] = W[i,j]
    
def rand_w(dist_name, dist_params, n, pz, upbound_ac, fn):
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
    ensure_upbound_ac(W, upbound_ac)
    print "saving matrix to " + fn
    mmwrite(fn, csr_matrix(W))

# main
if __name__ == '__main__':
    if len(argv) < 6:
        args = "<pz> <dist_name> <dist_params> <upbound_ac> <n1> [<n2> <n3>...]"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    pz = float(argv[1])
    dist_name = argv[2]
    dist_params = argv[3]
    upbound_ac = float(argv[4])
    ns  = map(int, argv[5:])
    fns = map(lambda n: str(n) + ".mtx", ns)
    for n,fn in zip(ns, fns):
        rand_w(dist_name, dist_params, n, pz, upbound_ac, fn)
