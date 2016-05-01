#!/usr/bin/env python

from sys import argv, exit
from test_util import get_lap, cond_err, parse_bool, eprint
    
# main
if __name__ == '__main__':
    if len(argv) < 3:
        args = "<is_lap> <file1> [<file2> ... ]"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    is_lap = parse_bool(argv[1])
    fns = argv[2:]
    for fn in fns:
        L = get_lap(fn, is_lap, fmt="dense")
        cn, re = cond_err(L)
        print("%-20s %.3E   %.3E" % (fn, cn, re))
