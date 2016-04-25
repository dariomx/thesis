########################################################
# Point to check: testing the algebraic connectivity
# (second smallest eigenvalue), is enough; no need to
# calculate all eigenvalues and check which ones are
# nearly zero.

# sample graph taken from
# https://en.wikipedia.org/wiki/Laplacian_matrix#Example
#
# with some tweaked weights
########################################################3

import networkx as nx
import scipy.sparse as sps
import scipy.linalg as spla
import unittest2 as ut
import sys
        
class TestDisconn1(ut.TestCase):
    def create_fully_conn_graph(self):
        G = nx.Graph()
        G.add_nodes_from(range(1,7))
        G.add_edge(1, 5, weight=1)
        G.add_edge(1, 2, weight=1)
        G.add_edge(5, 2, weight=1)
        G.add_edge(2, 3, weight=1)
        G.add_edge(5, 4, weight=1)
        G.add_edge(4, 3, weight=1)
        G.add_edge(6, 4, weight=1)
        return G;

    def create_fully_disconn_graph2(self):
        G = nx.Graph()
        G.add_nodes_from(range(1,7))
        G.add_edge(1, 5, weight=1)
        G.add_edge(1, 2, weight=1)
        G.add_edge(5, 2, weight=1)
        G.add_edge(2, 3, weight=1)
        G.add_edge(5, 4, weight=0)
        G.add_edge(4, 3, weight=0)
        G.add_edge(6, 4, weight=1)
        return G;

    def create_fully_disconn_graph3(self):
        G = nx.Graph()
        G.add_nodes_from(range(1,7))
        G.add_edge(1, 5, weight=0)
        G.add_edge(1, 2, weight=0)
        G.add_edge(5, 2, weight=1)
        G.add_edge(2, 3, weight=1)
        G.add_edge(5, 4, weight=0)
        G.add_edge(4, 3, weight=0)
        G.add_edge(6, 4, weight=1)
        return G;
    
    def test_full_conn_laplacian(self):
        L = sps.lil_matrix((6, 6))
        L.setdiag([2, 3, 2, 3, 3, 1])
        minus_ones = [(0,1), (0,4),
                      (1,0), (1,2),(1,4),
                      (2,1), (2,3),
                      (3,2), (3,4), (3,4), (3,5),
                      (4,0), (4,1), (4,3),
                      (5,3)]
        for (i,j) in minus_ones:
            L[i, j] = -1
        G = self.create_fully_conn_graph()
        diff = L - nx.laplacian_matrix(G)
        deb("dead1 =\n%s", L.toarray())
        deb("dead2 = \n%s\n", nx.laplacian_matrix(G).toarray())
        deb("dead3 = %s\n%s\n", diff.nonzero(), diff.toarray())        
        self.assertTrue(sum(map(len, diff.nonzero())) == 0)

    def test_algconn_fully_conn(self):
        G = self.create_fully_conn_graph()
        ac = nx.algebraic_connectivity(G)
        self.assertTrue(ac > 0)
        l, Q = spla.eigh(nx.laplacian_matrix(G).todense())
        deb("ac = %.16f\nl0 = %.16f\nl =\n%s\nQ =\n%s\n",
            ac, l[1], l, Q)
        self.assertTrue(is_nearly_zero(l[0]))
        self.assertTrue(is_nearly_zero(ac - l[1]))

    def test_algconn_fully_disconn2(self):
        G = self.create_fully_disconn_graph2()
        ac = nx.algebraic_connectivity(G)
        l, Q = spla.eigh(nx.laplacian_matrix(G).todense())
        deb("ac = %.16f\nl =\n%s\nQ =\n%s\n", ac, l, Q)
        self.assertFalse(ac > 0)
        self.assertEqual(len(filter(is_nearly_zero, l)), 2)
        self.assertTrue(is_nearly_zero(l[0]))
        self.assertTrue(is_nearly_zero(l[1]))

    def test_algconn_fully_disconn3(self):
        G = self.create_fully_disconn_graph3()
        ac = nx.algebraic_connectivity(G)
        l, Q = spla.eigh(nx.laplacian_matrix(G).todense())
        deb("ac = %+.16f\nl0 = %+.16f\n" +
            "df = %+54.53f\nzt = %+54.53f\n" +
            "l =\n%s\nQ =\n%s\n",
            ac, l[1], ac - l[1], zero_tol, l, Q)
        self.assertFalse(ac > 0)        
        self.assertTrue(is_nearly_zero(ac - l[1]))        
        self.assertEqual(len(filter(is_nearly_zero, l)), 3)
        self.assertTrue(is_nearly_zero(l[0]))
        self.assertTrue(is_nearly_zero(l[1]))
        
# main
debug = False
#debug = True
eps = 2.220446049250313E-16
zero_tol = 100000.0 * eps

def log(fmt, *args):
    sys.stdout.write(fmt % args)
    sys.stdout.flush()

def deb(fmt, *args):
    global debug
    if debug:
        log(fmt, *args)

def is_nearly_zero(x):
    return True if abs(x) < zero_tol else False

if __name__ == '__main__':
    suite = ut.TestSuite()
    suite.addTest(TestDisconn1('test_algconn_fully_conn')) 
    suite.addTest(TestDisconn1('test_algconn_fully_disconn2'))
    suite.addTest(TestDisconn1('test_algconn_fully_disconn3')) 
    ut.TextTestRunner().run(suite)
