#!/usr/bin/env python

from sys import argv, exit
from itertools import izip
import matplotlib.pylab as plt
import numpy as np
from test_util import eprint

# main
if __name__ == '__main__':
    n = len(argv)
    if n < 3:
        eprint("Usage: %s <name1> <file1> [<name2> <file2>]\n")
        exit(1)
    names = [argv[i] for i in xrange(1,n,2)]
    files = [argv[i] for i in xrange(2,n,2)]
    for name,file in izip(names, files):
        y = np.fromfile(file, sep=" ")
        plt.plot(y, label=name)
    plt.legend(loc='upper left')
    plt.show()

