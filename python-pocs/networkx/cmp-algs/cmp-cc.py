#!/usr/bin/env python

from sys import argv, exit
from fiedler import fiedler_vector
from test_util import eprint, get_lap
from test_util import take_time
from scipy.sparse.csgraph import connected_components as cc
from scipy.linalg import eigh
import numpy as np

# calculates application labels
# (splitting in 2 groups per snd zero eigenvector)
def calc_app_lab(fv):
    fmin, fmax = np.min(fv), np.max(fv)
    fvlimit = fmin + (fmax-fmin)/1000.0
    indic = lambda x: 0 if x < fvlimit else 1
    vindic = np.vectorize(indic)
    return vindic(fv)

def print_nzgrp(pfx, lab):
    n = len(lab)
    nzg = [x[0] for x in zip(xrange(n), lab) if x[1]>0]
    eprint("%s non zero group = %s" % (pfx,nzg))
    
def validate_cc(L, fv):
    (ncc, cclab), time = take_time(lambda: cc(L))
    eprint("ncc = %d, took %10.8f" % (ncc,time))
    print_nzgrp("cc", cclab)
    print_nzgrp("app", calc_app_lab(fv))
    
# main
if __name__ == '__main__':
    if len(argv) != 2:
        args = "<mat_file>"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    L = get_lap(argv[1], "dense")
    _, V = eigh(L, eigvals=(1,3), overwrite_a=True)
    validate_cc(L, V[:,1])
