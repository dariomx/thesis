#!/usr/bin/env python

from sys import argv
import numpy as np
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
from scipy.linalg import orth
import scipy.stats
from test_util import eprint

def get_rand_dist(dist_name, dist_params):
    dist = getattr(scipy.stats, dist_name)
    dps = map(float, dist_params.split(","))
    return lambda s: dist.rvs(*dps[:-2], loc=dps[-2], scale=dps[-1], size=s)

# returns random matrix with given eigenvalues
# assumes that rand_dist produces non-singular matrices
def rand_mat_eigv(rand_dist, n, eigv, fn):
    print("random generator will use dist %s" % dist_name)
    D = np.diag(eigv)
    Q = orth(rand_dist((n,n)))
    M = np.dot(Q, np.dot(D, Q.T))
    print "saving matrix to " + fn
    mmwrite(fn, csr_matrix(M))

# get the desired eigenvalues: heuristic is to get something similar
# to the laplacian of a graph
# smallest eigenval = 0
# second smallest = l2
# rest = small variations from l2 (use given random dist generator)
def get_eigvals(rand_dist, n, l2):
    eigv = np.sort(rand_dist(n))
    eigv[0:1] = 0
    eigv[1:] += l2
    return eigv
    
# main
if __name__ == '__main__':
    if len(argv) < 4:
        args = "<dist_name> <dist_params> l2 <n1> [<n2> <n3>...]"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    dist_name = argv[1]
    dist_params = argv[2]
    l2 = float(argv[3])
    ns  = map(int, argv[4:])
    fns = map(lambda n: str(n) + ".mtx", ns)
    rand_dist = get_rand_dist(dist_name, dist_params)
    for n,fn in zip(ns, fns):
        eigv = get_eigvals(rand_dist, n, l2)
        rand_mat_eigv(rand_dist, n, eigv, fn)
