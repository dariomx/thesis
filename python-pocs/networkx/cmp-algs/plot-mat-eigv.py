#!/usr/bin/env python

from sys import argv, exit
from itertools import izip
import matplotlib.pylab as plt
import numpy as np
from test_util import eprint, get_lap, get_eigvals, is_lap
from test_util import save_plot, get_plot_axis

# main
if __name__ == '__main__':
    n = len(argv)
    if n < 4:
        eprint("Usage: %s <img_file> <name1> <file1> [<name2> <file2>]\n")
        exit(1)
    img_file = argv[1]
    names = [argv[i] for i in xrange(2,n,2)]
    files = [argv[i] for i in xrange(3,n,2)]
    ax = get_plot_axis()
    for name, mat_file in izip(names, files):
        L = get_lap(mat_file, is_lap(mat_file), fmt="dense")
        eigv = get_eigvals(L)
        eprint("plotting eigenvalues of %s" % mat_file)
        ax.plot(eigv, label=name)
    ax.legend(loc='upper left')
    eprint("saving plots into %s" % img_file)
    save_plot(ax, "eigenvalues", img_file)
