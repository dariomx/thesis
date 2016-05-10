#!/usr/bin/env python

from sys import argv, exit
from test_util import load_weights, eprint
from test_util import save_plot, get_plot_axis, get_pdf_dist

# main
if __name__ == '__main__':
    if len(argv) not in (4,6):
        args = "<bins> <mat_file> <img_file> [<dist_name> <dist_params>]"
        eprint (("\nUsage: %s " + args + "\n") % argv[0])
        exit(1)
    bins = int(argv[1])
    mat_file = argv[2]
    img_file = argv[3]
    dist_name = argv[4] if len(argv)==6 else None
    dist_params = argv[5] if len(argv)==6 else None
    _, ws = load_weights(mat_file)
    eprint("building histogram ... ")
    ax = get_plot_axis()
    _, xs, _ = ax.hist(ws, bins=bins, color="b", normed=True)
    if dist_name is not None:
        eprint("plotting pdf of distribution %s ... " % dist_name)
        pdf = get_pdf_dist(dist_name, dist_params)
        ax.plot(xs, pdf(xs), label=dist_name, color="r")
    eprint("saving histogram to %s ... " % img_file)
    title = "histogram of matrix %s" % mat_file
    save_plot(ax, title, img_file)
