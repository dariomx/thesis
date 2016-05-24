#!/usr/bin/env python

from sys import argv, exit, exc_info
from test_util import get_lap, parse_bool, eprint, load_weights
from test_util import save_plot, get_plot_axis
import matplotlib.pyplot as plt
import scipy
import scipy.stats
import numpy as np

def fit_data(xs, ys, ws, dist_name):
    eprint("trying to fit data against dist %s" % dist_name)
    dist = getattr(scipy.stats, dist_name)
    try:
        param = dist.fit(ws)
        ys_fit = dist.pdf(xs, *param[:-2], loc=param[-2], scale=param[-1])
        fit_err = np.linalg.norm(ys - ys_fit[0:-1])
        eprint("%s gave fit_err of %f" % (dist_name, fit_err))
        return dist_name, param, ys_fit, fit_err
    except:
        args = (dist_name, str(exc_info()[0]))
        eprint("fitting of data against %s failed: %s" % args)
        return dist_name, None

def plot_best_fits(xs, fits, top_n, colors):
    is_valid = lambda x: x[1] is not None and not(np.isnan(x[3]))
    valid_fits = filter(is_valid, fits)
    best_fits = zip(xrange(top_n), 
                    sorted(valid_fits, key=lambda x: x[3])[0:top_n])
    for i, (dist_name, param, ys_fit, fit_err) in best_fits:
        args = (dist_name, top_n, fit_err)
        eprint("%s is amont top %d with fit_err=%s" % args)
        args = (dist_name, param)
        eprint("%s params = %s" % args)
        plt.plot(xs, ys_fit, label=dist_name, color=colors[i])
    plt.xlim(0,1)

def get_file_column(i, fn):
    return [l.split()[i] for l in open(fn, "r").read().splitlines()]

# main
if __name__ == '__main__':
    if len(argv) < 6:
        args = "<mat_file> <bins> <dist_names_file> <topn> <colors_file> <img_file>"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    mat_file = argv[1]
    bins = int(argv[2])
    dist_names = get_file_column(0, argv[3])
    topn = int(argv[4])
    colors = get_file_column(0, argv[5])
    img_file = argv[6]
    eprint("will try to fit against %d dists " % len(dist_names))
    _, _, n, ws = load_weights(mat_file)
    eprint("building histogram ... ")
    ax = get_plot_axis()
    ys, xs, _ = ax.hist(ws, bins=bins, color='blue', normed=True)
    fitd = lambda dn: fit_data(xs, ys, ws, dn)
    plot_best_fits(xs, map(fitd, dist_names), topn, colors)
    ax.legend(loc='upper right')
    eprint("saving plots to %s ..." % img_file)
    save_plot(ax, "best pdf fit for %s" % mat_file, img_file)
