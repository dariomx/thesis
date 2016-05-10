#!/usr/bin/env python

from sys import argv
from test_util import load_graph
import networkx as nx

# draws a graph from its laplacian
def draw(L_file, img_file):
    G  = load_graph(L_file, is_lap=True, zeros=True)
    aG = nx.to_agraph(G)
    aG.layout(prog="sfdp")
    aG.draw(img_file)
    
# main
if __name__ == '__main__':
    draw(argv[1], argv[2])
    
    
