#!/usr/bin/env python

from sys import argv, exit
from test_util import eprint, get_lap, conv_mat
from test_util import take_time
from scipy.sparse.csgraph import connected_components as cc
from scipy.io import mmwrite
from scipy.sparse import csc_matrix, coo_matrix
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
    L = get_lap(mat_file, "dense")
    eprint("Recovering original weights ...")
    W = -L
    np.fill_diagonal(W, 0)
    W = conv_mat(W, "csr")
    (ncc, cclab), time = take_time(lambda: cc(W, directed=False))
    eprint("ncc = %d, took %10.8f" % (ncc,time))
    if ncc > 1 and out_dir is not None:
        W = conv_mat(W, "dense")
        ks = get_nzgrp_idxs(cclab)
        eprint("removing non-zero-group nodes %s ..." % str(ks))
        W = np.delete(W, ks, axis=0)
        W = np.delete(W, ks, axis=1)
        eprint("checking for symmetry ...")
        assert (W.transpose() == W).all()
        eprint("recomputing laplacian ...")
        L = -W
        d = np.abs(np.sum(W, axis=0))
        L[np.diag_indices(L.shape[0])] = d
        new_mat_file = out_dir + "/" + str(L.shape[0])
        new_mat_file += "-" + basename(mat_file).split(".")[0]
        new_mat_file += "" if new_mat_file.endswith("-lap") else "-lap"
        new_mat_file += ".mtx"
        args = (str(L.shape), new_mat_file)
        eprint("converting back to sparse format (coo) ... ")
        L = coo_matrix(L)
        eprint("saving new laplacian of %s -> %s" % args)        
        mmwrite(new_mat_file, L)
