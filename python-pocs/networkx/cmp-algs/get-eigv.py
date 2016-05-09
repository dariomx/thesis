#!/usr/bin/env python

from sys import argv, exit
from scipy.linalg import eigh
from test_util import parse_bool, eprint, get_lap, get_eigvals
from numpy import savetxt

# main
if __name__ == '__main__':
    if len(argv) != 4:
        args = "<is_lap> <mat_file> <eigv_file>"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    is_lap = parse_bool(argv[1])
    mat_file = argv[2]
    eigv_file = argv[3]
    L = get_lap(mat_file, is_lap, fmt="dense")
    eigv = get_eigvals(L)
    savetxt(eigv_file, eigv)
