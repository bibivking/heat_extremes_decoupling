#!/usr/bin/env python

"""
Plot a histogram showing the distribution of positive/negative LE slopes

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (20.04.2018)"
__email__ = "mdekauwe@gmail.com"

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import seaborn as sns
import sys

def main(fname):

    df = pd.read_csv(fname)

    """
    p_events = df[df.slope > 0.0].doy
    p_events = np.floor(p_events)
    unique, counts = np.unique(p_events, return_counts=True)
    unique = [i.astype(int) for i in unique]
    df1 = pd.DataFrame(unique, columns = ['doy'])
    df1['counts'] = counts

    n_events = df[df.slope < 0.0].doy
    n_events = np.floor(n_events)

    unique, counts = np.unique(n_events, return_counts=True)
    unique = [i.astype(int) for i in unique]
    df2 = pd.DataFrame(unique, columns = ['doy'])
    df2['counts'] = counts

    sns.distplot(df1.doy, bins=30, kde=False, hist=True)
    sns.distplot(df2.doy, bins=30, kde=False, hist=True)
    plt.show()
    """
    p_events = df[df.slope > 0.0].doy
    p_events = np.floor(p_events)

    n_events = df[df.slope < 0.0].doy
    n_events = np.floor(n_events)
    unique, counts = np.unique(n_events, return_counts=True)

    
    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    ax1 = fig.add_subplot(111)
    ax1.hist(p_events, bins=30, density=False, ls='-', lw=1, facecolor="None",
             edgecolor="black", label="Positive slope")
    ax1.hist(n_events, bins=30, color="green", density=False,
             label="Negative slope", alpha=0.4)

    ax1.legend(numpoints=1, loc="best", frameon=False)
    ax1.set_ylabel("Count")
    ax1.set_xlabel("Day of year")

    fdir = "/Users/mdekauwe/Dropbox/fluxnet_heatwaves_paper/figures/figs"
    ofname = "timing_Qle.pdf"
    fig.savefig(os.path.join(fdir, ofname),
                bbox_inches='tight', pad_inches=0.1)

    #plt.show()


if __name__ == "__main__":

    data_dir = "data"
    fname = os.path.join(data_dir, "ozflux_slopes_Qle.csv")

    main(fname)
