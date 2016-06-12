#!/usr/bin/env python

from sys import argv, exit
from datetime import datetime
from scipy.linalg import eigh
from fiedler import fiedler_vector
from test_util import relres, parse_bool, eprint, get_lap, conv_mat
from test_util import get_ac_upbound, is_nearly_zero, invsign
from test_util import take_time
from fiedler_power import fiedler_pm, fiedler_invpow, fiedler_ship
from fiedler_power import fiedler_suip, fiedler_surqi, fiedler_rsuip
from fiedler_power import fiedler_suipc
from numpy import sign

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
    elif method == "suipc":
        return fiedler_suipc(L)    
    elif method == "rsuip":
        return fiedler_rsuip(L)    
    elif method == "surqi":
        return fiedler_surqi(L)    
    else: 
        return fiedler_vector(L, method=method)

def test_fiedler(L, method):
    try:
        f = lambda: calc_fiedler(L, method) 
        (ac, fv), time = take_time(f)
        res = relres(L, ac, fv)
        return ac, fv, time, res        
    except:
        eprint("Error while calculating %s", method)
        raise
        return None

def relname(fn):
    return fn.split("/")[-1]

def validate_fv(ac, fv, gac, gfv):
    ac_ok = (sign(ac) == sign(gac) and
             is_nearly_zero(ac) ==
             is_nearly_zero(gac))
    fv_signok = 0
    for x,y in zip(invsign(fv, gfv), gfv):
        if sign(x) != sign(y):
            #eprint("fv diff: %f %f" % (x,y))
            None
        else:
            fv_signok += 1
    fv_ok = fv_signok / float(fv.shape[0])
    args = (ac_ok, fv_signok, fv_ok)
    return "ac_ok=%s fv_signok=%d,fv_ok=%.7f" % args

def validate_cc(L, fv):
    ncc, cclab = cc(L)
    eprint("ncc = %d" % ncc)        
    

def test_methods(methods, fmt, fns):
    gac, gfv, valstr = None, None, ""
    for fn in fns:
        L = get_lap(fn, fmt)
        for met in methods:
            ret = test_fiedler(L, met)
            if ret is None:
                continue
            ac, fv, time, res = ret
            if met == "mr3":
                gac, gfv = ac, fv
                valstr = ""
            elif gac is not None and gfv is not None:
                valstr = validate_fv(ac, fv, gac, gfv)
            args = (fn, met, time, ac, res, valstr)
            if not (met == "mr3" and verif):
                print("%-30s %-10s %10.8f  %.16f  %.3E %s" % args)
    
# main
if __name__ == '__main__':
    if len(argv) < 5:
        args = "<method> <fmt> <verif> <file1> [<file2> ... ]"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    method_list = argv[1]
    fmt = argv[2]
    verif = parse_bool(argv[3])
    fns = argv[4:]
    if method_list == "all":
        methods = ["mr3", "lanczos_susi", "suip"]
    else:
        methods = method_list.split(",")
    if verif and "mr3" not in methods:
        methods = ["mr3"] + methods
    test_methods(methods, fmt, fns)
            
