#!/usr/bin/env python

from sys import argv
from scipy.io import mmwrite
from test_util import eprint, load_mat

# converts a (presumably symmetric) matrix into mm general format 
def symm2gen(in_file, out_file):
    if in_file != out_file:
        print "converting %s -> %s" % (in_file, out_file)
        sparse = load_mat(in_file)
        mmwrite(target=out_file,
                a=sparse,
                comment="miau",
                field="real",
                symmetry="general")
        print "converted %s -> %s" % (in_file, out_file)
    else:
        print "ignoring conversion, gave same files"
        

# main
if __name__ == '__main__':
    if len(argv) != 3:
        args = "<in_file> <out_file>"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    else:
        symm2gen(argv[1], argv[2])
    
