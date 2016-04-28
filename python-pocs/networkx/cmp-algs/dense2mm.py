#!/usr/bin/env python

from sys import argv
from numpy import loadtxt
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from test_util import parse_bool

# converts a file in dense (lil) format into a new file in matrix-market
# format (sparse flavor, with an optional symmetric flag)
def dense2mm(dense_file, sparse_file, symmetric):
    if dense_file != sparse_file:
        print "converting %s -> %s" % (dense_file, sparse_file)
        dense = loadtxt(dense_file)
        sparse = csr_matrix(dense)
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
    dense2mm(argv[1], argv[2], parse_bool(argv[3]))
    
