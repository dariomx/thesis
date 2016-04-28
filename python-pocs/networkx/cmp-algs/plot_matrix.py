#!/usr/bin/env python

from sys import argv
import matplotlib.pylab as plt
from test_util import load_mat

# main
if __name__ == '__main__':
    M_file = argv[1]
    M = load_mat(M_file, False)
    plt.spy(M)
    plt.show()

