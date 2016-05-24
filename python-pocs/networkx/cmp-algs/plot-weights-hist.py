#!/usr/bin/env python

from sys import argv, exit
from test_util import load_weights, eprint
from test_util import save_plot, get_plot_axis, get_pdf_dist
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plot_2d_hist(ws, ax, bins, dist_name, dist_params):
    _, xs, _ = ax.hist(ws, bins=bins, color="b", normed=True)
    if dist_name is not None:
        eprint("plotting pdf of distribution %s ... " % dist_name)
        pdf = get_pdf_dist(dist_name, dist_params)
        ax.plot(xs, pdf(xs), label=dist_name, color="r")    

def plot_3d_hist(ix, iy, ws, ax):
    n = len(ws)
    xpos = ix
    ypos = iy
    zpos = np.zeros(n)
    dx = np.ones(n) * 0.50
    dy = np.ones(n) * 0.50
    dz = ws
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='blue', zsort='average')
    
# main
if __name__ == '__main__':
    if len(argv) not in (5,7) or argv[4] not in ("2d","3d"):
        args = "<bins> <mat_file> <img_file> <2d|3d> [<dist_name> <dist_params>]"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    bins = int(argv[1])
    mat_file = argv[2]
    img_file = argv[3]
    hist_type = argv[4]
    dist_name = argv[5] if len(argv)==7 else None
    dist_params = argv[6] if len(argv)==7 else None
    ix, iy, _, ws = load_weights(mat_file)
    eprint("building %s histogram ... " % hist_type)
    ax = get_plot_axis(hist_type == "3d")    
    if hist_type == "2d":
        plot_2d_hist(ws, ax, bins, dist_name, dist_params)
    else:
        plot_3d_hist(ix, iy, ws, ax)
    eprint("saving histogram to %s ... " % img_file)
    title = "histogram of matrix %s" % mat_file
    save_plot(ax, title, img_file)
