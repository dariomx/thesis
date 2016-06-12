#!/usr/bin/env python

from sys import argv, exit
from test_util import eprint, get_lap, conv_mat
from test_util import take_time
from scipy.sparse.csgraph import connected_components as cc
from scipy.io import mmwrite
from scipy.sparse import csc_matrix
import numpy as np
from os.path import basename

def get_nzgrp_idxs(lab):
    n = len(lab)
    return [x[0] for x in zip(xrange(n), lab) if x[1]>0]

if __name__ == '__main__':
    if len(argv) < 2:
        args = "<mat_file> [<out_dir>]"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    mat_file = argv[1]
    out_dir = None if len(argv)==2 else argv[2]
    L = get_lap(mat_file, "csc")
    (ncc, cclab), time = take_time(lambda: cc(L))
    eprint("ncc = %d, took %10.8f" % (ncc,time))
    if out_dir is not None:
        L = conv_mat(L, "dense")
        for k in get_nzgrp_idxs(cclab):
            eprint("removing non-zero-group node %d" % k)
            np.fill_diagonal(L, L.diagonal() - L[:,k])
            L = np.delete(L, k, axis=0)
            L = np.delete(L, k, axis=1)
        L = csc_matrix(L)    
        new_mat_file = out_dir + "/" + str(L.shape[0])
        new_mat_file += "-" + basename(mat_file)
        args = (str(L.shape), new_mat_file)
        eprint("saving new laplacian of %s -> %s" % args)
        mmwrite(new_mat_file, L)
