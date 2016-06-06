#!/usr/bin/env python

from test_util import load_mat, eprint
from scipy.linalg import solve
from scipy.linalg import norm
from scipy.io import mmwrite
from numpy import dot
from sys import argv

# main
if __name__ == '__main__':
    if len(argv) < 5:
        args = "<A_file> <b_file> <x_file> <x_file_calc>"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    A = load_mat(argv[1]).toarray()
    b = load_mat(argv[2]).toarray()
    x = load_mat(argv[3]).toarray()
    res = norm(dot(A, b) - x) / norm(b)
    print("in res = %.14f" % res)
    x = solve(A, b)
    res = norm(dot(A, x) - b) / norm(b)
    print("calc res = %.14f" % res)
    mmwrite(argv[4], x);
    print("calc x saved to %s" % argv[4])
