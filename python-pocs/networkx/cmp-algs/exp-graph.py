#!/usr/bin/env python

from sys import argv
from os.path import isfile
from test_util import load_graph, parse_bool, is_lap
import networkx as nx

# exports a graph from its weight or laplacian matrix into desired format.
def exp_graph(M_file, is_lap, zeros, G_file, fmt):
    assert not isfile(G_file), "Graph output file already exist"
    G = load_graph(M_file, is_lap, zeros)
    print "writing graph -> %s" % G_file
    if fmt == "gexf":
        nx.write_gexf(G, G_file, prettyprint=True)
    elif fmt == "dot":
        nx.write_dot(G, G_file)
    else:
        raise ValueError("Unsupported format " + fmt)
    
# main
if __name__ == '__main__':
    M_file = argv[1]
    zeros  = parse_bool(argv[2])
    G_file = argv[3]
    fmt    = argv[4]
    exp_graph(M_file, is_lap(M_file), zeros, G_file, fmt)
