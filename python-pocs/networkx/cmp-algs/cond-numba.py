#!/usr/bin/env python

from sys import argv, exit
from test_util import get_lap, cond_err, parse_bool, eprint, is_lap
    
# main
if __name__ == '__main__':
    if len(argv) < 2:
        args = "<file1> [<file2> ... ]"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    fns = argv[1:]
    for fn in fns:
        L = get_lap(fn, is_lap(fn), fmt="dense")
        cn, re = cond_err(L)
        print("%-20s %.3E   %.3E" % (fn, cn, re))
