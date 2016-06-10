#!/usr/bin/env python

from sys import argv, exit
from test_util import eprint, get_lap
from scipy.io import mmwrite
from os.path import basename

# main
if __name__ == '__main__':
    if len(argv) < 3:
        args = "<eigv_outdir> <mat_infile>*"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    outdir = argv[1]
    for mat_file in argv[2:]:
        eprint("computing laplacian of %s" % mat_file) 
        L = get_lap(mat_file, fmt="csc")
        prefix = basename(mat_file).split(".")[0]
        lap_file = outdir + "/" + prefix + "-lap.mtx"
        mmwrite(lap_file, L)
