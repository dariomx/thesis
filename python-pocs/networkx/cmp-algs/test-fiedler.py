#!/usr/bin/env python

from sys import argv, exit
from datetime import datetime
from scipy.linalg import eigh
from fiedler import fiedler_vector
from test_util import relres, parse_bool, eprint, get_lap, conv_mat
from test_util import get_ac_upbound
from fiedler_power import fiedler_pm, fiedler_invpow, fiedler_ship
from fiedler_power import fiedler_suip, fiedler_surqi 


def calc_fiedler(L, method):
    if method == "mr3":
        Ld = conv_mat(L, "dense")
        ls, vs = eigh(Ld, eigvals=(1,1), overwrite_a=True)
        ac = ls[0]
        fv = vs[:,0]
        return ac, fv
    elif method == "pm":
        return fiedler_pm(L)
    elif method == "invpow":
        return fiedler_invpow(L)
    elif method == "ship":
        return fiedler_ship(L)    
    elif method == "suip":
        return fiedler_suip(L)
    elif method == "surqi":
        return fiedler_surqi(L)    
    else: 
        return fiedler_vector(L, method=method)

def test_fiedler(L, method):
    start = datetime.now()
    ac, fv = calc_fiedler(L, method)
    end = datetime.now()
    time = (end - start).total_seconds()
    res = relres(L, ac, fv)
    return ac, fv, time, res

def relname(fn):
    return fn.split("/")[-1]

# main
if __name__ == '__main__':
    if len(argv) < 4:
        args = "<method> <fmt> <file1> [<file2> ... ]"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    method_list = argv[1]
    fmt = argv[2]
    fns = argv[3:]
    if method_list == "all":
        methods = ["mr3", "lanczos", "lobpcg", "tracemin_pcg"]
    else:
        methods = method_list.split(",")
    for fn in fns:
        L = get_lap(fn, fmt)
        for met in methods:
            ac, fv, time, res = test_fiedler(L, met)
            args = (fn, met, time, ac, res)
            print("%-30s %-10s %10.8f  %.3E  %.3E" % args)
            #eprint("%s" % (fv))
            
