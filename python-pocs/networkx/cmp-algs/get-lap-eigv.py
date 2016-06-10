#!/usr/bin/env python

from sys import argv, exit
from scipy.linalg import eigh
from test_util import eprint, get_lap, get_eigvals, is_nearly_zero, take_time
from numpy import savetxt
from os.path import basename

# main
if __name__ == '__main__':
    if len(argv) < 3:
        args = "<eigv_outdir> <mat_infile>*"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    outdir = argv[1]
    discon = 0
    for mat_file in argv[2:]:
        eprint("computing eigenvalues for %s" % mat_file)        
        L = get_lap(mat_file, fmt="dense")
        eigv, time = take_time(lambda: get_eigvals(L))
        eprint("computing all eigenvalues took %10.8f secs" % time)
        ac = eigv[1]
        if is_nearly_zero(ac):
            discon += 1
            args = (mat_file, discon)
            eprint("graph for %s is disconnected (seen %d so far)" % args)
        else:
            eprint("ac = %.16f" % ac)
        eigv_file = outdir + "/eigv." + basename(mat_file).split(".")[0]
        eprint("saving eigenvalues to %s" % eigv_file)
        savetxt(eigv_file, eigv)
    eprint("found %d disconnected graphs" % discon)
