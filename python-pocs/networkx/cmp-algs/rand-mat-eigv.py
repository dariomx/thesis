#!/usr/bin/env python

from sys import argv
import numpy as np
from test_util import eprint, get_rand_dist, rand_mat_eigv
from test_util import get_lap, get_eigvals
from scipy.io import mmwrite

# main
if __name__ == '__main__':
    if len(argv) != 5:
        args = "<dist_name> <dist_params> <n> <mat_file>"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    dist_name = argv[1]
    dist_params = argv[2]
    n  = int(argv[3])
    mat_file = argv[4]
    eprint("getting eigenvalues from matrix %s ..." % mat_file)
    L = get_lap(mat_file, fmt="dense")
    eigv = get_eigvals(L)
    eprint("setting eigenvalues on random matrix (dist %s) ... " % dist_name)
    rand_dist = get_rand_dist(dist_name, dist_params)    
    out_file = str(n) + "-lap.mtx"
    M = rand_mat_eigv(rand_dist, n, eigv)
    eprint("saving matrix to " + out_file)
    mmwrite(out_file, M)

