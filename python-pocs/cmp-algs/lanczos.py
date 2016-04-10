#!/usr/bin/env python

from sys import argv
from fiedler import fiedler_vector
from util import test

def calc(L):
    return fiedler_vector(L, method="lanczos")

# main
if __name__ == '__main__':
    test(argv, calc)
    
