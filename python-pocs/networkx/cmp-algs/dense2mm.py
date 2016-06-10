#!/usr/bin/env python

from sys import argv
from numpy import loadtxt
from scipy.io import mmwrite
from scipy.sparse import csc_matrix
from test_util import parse_bool, eprint
from os.path import basename

# converts a file in dense (lil) format into a new file in matrix-market
# format (sparse flavor, with an optional symmetric flag)
def dense2mm(dense_file, sparse_file, symmetric):
    if dense_file != sparse_file:
        print "converting %s -> %s" % (dense_file, sparse_file)
        dense = loadtxt(dense_file)
        sparse = csc_matrix(dense)
        symmetry = "symmetric" if symmetric else "general"
        mmwrite(target=sparse_file,
                a=sparse,
                comment="miau",
                field="real",
                symmetry=symmetry)
    else:
        print "ignoring conversion, gave same files"
        

# main
if __name__ == '__main__':
    if len(argv) < 4:
        args = "<symmetric> <sparse_dir> <dense_file>*"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    symmetric = parse_bool(argv[1])
    sparse_dir = argv[2]
    for dense_file in argv[3:]:
        prefix = basename(dense_file).split(".")[0]
        sparse_file = sparse_dir + "/" + prefix + ".mtx"
        dense2mm(dense_file, sparse_file, symmetric)
    
