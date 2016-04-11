#!/usr/bin/env python

from sys import argv
from os.path import isfile
from test_util import load_graph
import networkx as nx

# exports a graph in gexf
def exp_gexf(L_file, G_file):
    assert not isfile(G_file), "Graph output file already exist"
    G = load_graph(L_file)
    print "writing graph -> %s" % G_file
    nx.write_gexf(G, G_file, prettyprint=True)
    
# main
if __name__ == '__main__':
    exp_gexf(argv[1], argv[2])
    
    
