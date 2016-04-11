#!/usr/bin/env python

from sys import argv
from numpy import loadtxt
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

# converts a file in dense (lil) format into a new file in matrix-market
# format (sparse flavor)
def dense2mm(dense_file, sparse_file):
    if dense_file != sparse_file:
        print "converting %s -> %s" % (dense_file, sparse_file)
        dense = loadtxt(dense_file)
        sparse = csr_matrix(dense)
        mmwrite(target=sparse_file,
                a=sparse,
                comment="laplacian of 867",
                field="real")
    else:
        print "ignoring conversion, gave same files"
        

# main
if __name__ == '__main__':
    dense2mm(argv[1], argv[2])
    
