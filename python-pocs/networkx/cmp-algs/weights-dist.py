#!/usr/bin/env python

from sys import argv, exit
from test_util import get_lap, get_weights, parse_bool, eprint
import matplotlib.pyplot as plt
import scipy
import scipy.stats
import numpy as np

def load_weights(is_lap, fn):
    L = get_lap(fn, is_lap, fmt="csr")
    eprint("extracting weights from %s ... " % fn)        
    ws = get_weights(L)
    nnz = np.squeeze(np.asarray(ws[ws > 0].T))
    args = (ws.shape[1], nnz.shape[0])
    eprint("extracted %d weights (%d nnz)" % args)
    return nnz.shape[0], nnz

def fit_data(xs, ys, ws, dist_name):
    eprint("trying to fit data against dist %s" % dist_name)
    dist = getattr(scipy.stats, dist_name)
    param = dist.fit(ws)
    ys_fit = dist.pdf(xs, *param[:-2], loc=param[-2], scale=param[-1])
    fit_err = np.linalg.norm(ys - ys_fit[0:-1])
    eprint("%s gave fit_err of %f" % (dist_name, fit_err))
    return dist_name, ys_fit, fit_err

def plot_best_fits(xs, fits, top_n):
    best_fits = sorted(fits, key=lambda x: x[2])[0:top_n]
    for dist_name, ys_fit, fit_err in best_fits:
        args = (dist_name, top_n, fit_err)
        eprint("%s is amont top %d with fit_err=%s" % args)
        plt.plot(xs, ys_fit, label=dist_name)
    plt.xlim(0,1)

def get_dist_names(fn):
    return [l.split()[0] for l in open(fn, "r").read().splitlines()]
    
# main
if __name__ == '__main__':
    if len(argv) < 4:
        args = "<is_lap> <file1> <bins> <dist_names_file>"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    is_lap = parse_bool(argv[1])
    fn = argv[2]
    bins = int(argv[3])
    dist_names = get_dist_names(argv[4])
    eprint("will try to fit against %d dists " % len(dist_names))
    n, ws = load_weights(is_lap, fn)
    eprint("building histogram ... ")
    ys, xs, _ = plt.hist(ws, bins=bins, color='w', normed=True)
    fitd = lambda dn: fit_data(xs, ys, ws, dn)
    plot_best_fits(xs, map(fitd, dist_names), 3)
    plt.legend(loc='upper right')
    plt.show()
