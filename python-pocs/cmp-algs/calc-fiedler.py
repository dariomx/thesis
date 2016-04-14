#!/usr/bin/env python

from sys import argv
from datetime import datetime
from scipy.linalg import eigh
from fiedler import fiedler_vector
from test_util import load_lap, load_mat
from test_util import cmp_ac_fv, lap, parse_bool

def calc_fiedler(L, method):
    if method == "mr3":
        ls, vs = eigh(L, eigvals=(1,1), overwrite_a=True)
        ac = ls[0]
        fv = vs[:,0]
        return ac, fv
    else: 
        return fiedler_vector(L, method=method)

def test(M_file, is_lap, method, dense, eigvals_file, fv_file):
    if is_lap:
        L, ac, fv = load_lap(M_file, dense, eigvals_file, fv_file)
    else:
        W = load_mat(M_file, dense)
        L = lap(W, dense)
        ac, fv = None, None
    start = datetime.now()
    cac, cfv = calc_fiedler(L, method)
    end = datetime.now()
    print "calc took %s" % (end - start)
    cmp_ac_fv(L, cac, cfv, ac, fv)
    
# main
if __name__ == '__main__':
    M_file = argv[1]
    is_lap = parse_bool(argv[2])
    method = argv[3]
    dense = parse_bool(argv[4]) if len(argv) >= 5 else False
    eigvals_file = argv[5] if len(argv) >=6 else None
    fv_file = argv[6] if len(argv) >= 7 else None
    test(M_file, is_lap, method, dense, eigvals_file, fv_file)
