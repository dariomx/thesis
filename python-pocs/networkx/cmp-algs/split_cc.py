#!/usr/bin/env python

from sys import argv, exit
from test_util import eprint, load_mat, lap, conv_mat, get_mat_cons
from test_util import take_time, parse_bool
from scipy.sparse.csgraph import connected_components as cc
from scipy.io import mmwrite
from scipy.sparse import csc_matrix, coo_matrix
import numpy as np
from os.path import basename
from itertools import izip

def get_cc(W):
    eprint("computing scc ...")
    ncc, cclab = cc(conv_mat(W,"csr"), directed=False)
    return ncc, cclab

def split_cc_dense(W, cclab):
    W = conv_mat(W, "dense")
    ks = np.nonzero(cclab)[0]
    eprint("removing non-zero-group nodes %s ..." % str(ks))
    W = np.delete(W, ks, axis=0)
    W = np.delete(W, ks, axis=1)
    return conv_mat(W,"csr")

def split_cc_sparse1(W, cclab):
    ks = np.nonzero(cclab)[0]
    eprint("removing non-zero-group nodes %s ..." % str(ks))
    def stays((d,i,j)):
        return i not in ks and j not in ks
    dd,ii,jj = izip(*filter(stays, izip(W.data, W.row, W.col)))
    W.data, W.row, W.col = np.array(dd), np.array(ii), np.array(jj)
    return W

# adapted from http://stackoverflow.com/questions/23966923/delete-columns-of-matrix-of-csr-format-in-python
def split_cc_sparse2(W, cclab):
    idx_del = np.nonzero(cclab)[0]
    eprint("removing non-zero-group nodes %s ..." % str(idx_del))
    keep_row = np.logical_not(np.in1d(W.row, idx_del))
    keep_col = np.logical_not(np.in1d(W.col, idx_del))
    keep = np.logical_and(keep_row, keep_col)
    W.data = W.data[keep]
    W.row = W.row[keep]
    W.col = W.col[keep]
    W.row -= np.less.outer(idx_del, W.row).sum(0)
    W.col -= np.less.outer(idx_del, W.col).sum(0)
    k = len(idx_del)
    W._shape = (W.shape[0] - k, W.shape[1] - k)
    return W

def save_lap(L, mat_file, out_dir):
    new_mat_file = out_dir + "/" + str(L.shape[0])
    new_mat_file += "-" + basename(mat_file).split(".")[0]
    new_mat_file += "" if new_mat_file.endswith("-lap") else "-lap"
    new_mat_file += ".mtx"
    args = (str(L.shape), new_mat_file)        
    eprint("saving new laplacian of %s -> %s" % args)        
    mmwrite(new_mat_file, L)    
    
if __name__ == '__main__':
    if len(argv) < 5:
        args = "<split> <navg> <save> <out_dir> <mat_file>*"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    split = argv[1]
    navg = int(argv[2])
    save = parse_bool(argv[3])
    out_dir = argv[4]
    split_cc = split_cc_sparse2 if split=="sparse" else split_cc_dense
    for mat_file in argv[5:]:
        W = load_mat(mat_file, "coo")
        (ncc, cclab), time_cc = take_time(lambda: get_cc(W), navg)
        if ncc > 1:
            W, time_split = take_time(lambda: split_cc(W, cclab), navg)
        else:
            time_split = 0
        eprint("recomputing laplacian ...")
        L, time_lap = take_time(lambda: lap(W, "csr"))
        args = (mat_file, ncc, time_cc, time_split, time_lap,
                time_cc + time_split + time_lap) 
        eprint("%s, %d, %10.8f, %10.8f, %10.8f, %10.8f" % args)
        if save:
            _, time_save = take_time(lambda: save_lap(L, mat_file, out_dir))
            eprint("saving matrix took %10.8f" % time_save)
