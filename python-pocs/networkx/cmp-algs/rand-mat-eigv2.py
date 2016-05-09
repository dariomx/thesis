#!/usr/bin/env python

from sys import argv
import numpy as np
from test_util import eprint, get_rand_dist, rand_mat_eigv

# main
if __name__ == '__main__':
    if len(argv) < 4:
        args = "<dist_name> <dist_params> <eigv_file>"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    dist_name = argv[1]
    dist_params = argv[2]
    eigv_file = argv[3]
    rand_dist = get_rand_dist(dist_name, dist_params)
    eigv = np.fromfile(eigv_file, sep=" ")
    n  = eigv.shape[0]    
    mat_file = str(n) + ".mtx"
    rand_mat_eigv(rand_dist, n, eigv, mat_file)
