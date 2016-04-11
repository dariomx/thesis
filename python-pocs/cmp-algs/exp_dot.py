#!/usr/bin/env python

from sys import argv
from os.path import isfile
from test_util import load_graph
import networkx as nx

# exports a graph in gexf
def exp_dot(L_file, G_file, zeros=True):
    assert not isfile(G_file), "Graph output file already exist"
    G = load_graph(L_file, zeros)
    print "writing graph -> %s" % G_file
    nx.write_dot(G, G_file)
    
# main
if __name__ == '__main__':
    exp_dot(argv[1], argv[2], argv[3])
