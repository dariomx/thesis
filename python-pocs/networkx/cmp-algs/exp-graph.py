#!/usr/bin/env python

from sys import argv
from os.path import isfile
from test_util import load_graph, parse_bool
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
    is_lap = parse_bool(argv[2])
    zeros  = parse_bool(argv[3])
    G_file = argv[4]
    fmt    = argv[5]
    exp_graph(M_file, is_lap, zeros, G_file, fmt)
