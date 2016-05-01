#!/usr/bin/env python

from sys import argv, exit
from datetime import datetime
from scipy.linalg import eigh
from fiedler import fiedler_vector
from test_util import relres, parse_bool, eprint, get_lap, conv_mat

def calc_fiedler(L, method):
    if method == "mr3":
        Ld = conv_mat(L, "dense")
        ls, vs = eigh(Ld, eigvals=(1,1), overwrite_a=True)
        ac = ls[0]
        fv = vs[:,0]
        return ac, fv
    else: 
        return fiedler_vector(L, method=method)

def test_fiedler(L, method):
    start = datetime.now()
    ac, fv = calc_fiedler(L, method)
    end = datetime.now()
    time = (end - start).total_seconds()
    res = relres(L, ac, fv)
    return ac, time, res
    
# main
if __name__ == '__main__':
    if len(argv) < 5:
        args = "<is_lap> <method> <fmt> <file1> [<file2> ... ]"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    is_lap = parse_bool(argv[1])
    method = argv[2]
    fmt = argv[3]
    fns = argv[4:]
    if method == "all":
        methods = ["mr3", "lanczos", "lobpcg", "tracemin_pcg"]
    else:
        methods = [method]
    for fn in fns:
        L = get_lap(fn, is_lap, fmt)
        for met in methods:
            ac, time, res = test_fiedler(L, met)
            print("%-20s %-15s %5.2f\t%-.3E\t%.3E" % (fn, met, time, ac, res))
