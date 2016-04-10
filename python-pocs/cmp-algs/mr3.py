#!/usr/bin/env python

from sys import argv
from scipy.linalg import eigh
from util import test

def calc(L):
    ls, vs = eigh(L, eigvals=(1,1), overwrite_a=True)
    ac = ls[0]
    fv = vs[:,0]
    return ac, fv

# main
if __name__ == '__main__':
    test(argv, calc)
