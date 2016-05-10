#!/usr/bin/env python

from sys import argv
import numpy as np
from scipy.io import mmwrite
from test_util import eprint, get_rand_dist, rand_mat_eigv

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
    fns = map(lambda n: str(n) + "-lap.mtx", ns)
    rand_dist = get_rand_dist(dist_name, dist_params)
    for n,fn in zip(ns, fns):
        eigv = get_eigvals(rand_dist, n, l2)
        M = rand_mat_eigv(rand_dist, n, eigv)
        print("saving matrix to " + fn)
        mmwrite(fn, M)
