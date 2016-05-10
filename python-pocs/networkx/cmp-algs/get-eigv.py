#!/usr/bin/env python

from sys import argv, exit
from scipy.linalg import eigh
from test_util import eprint, get_lap, get_eigvals, is_lap
from numpy import savetxt

# main
if __name__ == '__main__':
    if len(argv) != 3:
        args = "<mat_file> <eigv_file>"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    mat_file = argv[1]
    eigv_file = argv[2]
    L = get_lap(mat_file, is_lap(mat_file), fmt="dense")
    eigv = get_eigvals(L)
    savetxt(eigv_file, eigv)
