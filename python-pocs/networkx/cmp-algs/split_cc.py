#!/usr/bin/env python

from sys import argv, exit
from test_util import eprint, load_mat, get_lap, conv_mat
from test_util import take_time, parse_bool
from scipy.sparse.csgraph import connected_components as cc
from scipy.io import mmwrite
from scipy.sparse import csc_matrix, coo_matrix
import numpy as np
from os.path import basename

def get_cc(W):
    eprint("computing scc ...")
    ncc, cclab = cc(W, directed=False)
    return ncc, cclab

def get_nzgrp_idxs(lab):
    n = len(lab)
    return [x[0] for x in zip(xrange(n), lab) if x[1]>0]

def split_cc(W, ncc, cclab):
    W = conv_mat(W, "dense")
    ks = get_nzgrp_idxs(cclab)
    eprint("removing non-zero-group nodes %s ..." % str(ks))
    W = np.delete(W, ks, axis=0)
    W = np.delete(W, ks, axis=1)
    #eprint("checking for symmetry ...")
    #assert (W.transpose() == W).all()
    eprint("recomputing laplacian ...")
    L = -W
    d = np.abs(np.sum(W, axis=0))
    L[np.diag_indices(L.shape[0])] = d
    return conv_mat(L, "csr")

def save_lap(L, mat_file, out_dir):
    new_mat_file = out_dir + "/" + str(L.shape[0])
    new_mat_file += "-" + basename(mat_file).split(".")[0]
    new_mat_file += "" if new_mat_file.endswith("-lap") else "-lap"
    new_mat_file += ".mtx"
    args = (str(L.shape), new_mat_file)        
    eprint("converting back to sparse format (coo) ... ")
    L = coo_matrix(L)
    eprint("saving new laplacian of %s -> %s" % args)        
    mmwrite(new_mat_file, L)    
    
if __name__ == '__main__':
    if len(argv) < 4:
        args = "<navg> <save> <out_dir> <mat_file>*"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    navg = int(argv[1])
    save = parse_bool(argv[2])
    out_dir = argv[3]
    for mat_file in argv[4:]:
        W = load_mat(mat_file, "csr")
        (ncc, cclab), time_cc = take_time(lambda: get_cc(W), navg=navg)
        if ncc > 1:
            L, time_split = take_time(lambda: split_cc(W, ncc, cclab), navg=navg)
        else:
            time_split = 0
            L = get_lap(W)
        args = (mat_file, ncc, time_cc, time_split)  
        eprint("%s, %d, %10.8f, %10.8f" % args)
        if save:
            save_lap(L, mat_file, out_dir)
