#!/usr/bin/env python

"""
For each of the OzFlux sites plot a histogram of the positive/negative slopes
for GPP

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (20.04.2018)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import re
import string

def main(fname):

    df = pd.read_csv(fname)
    sites = np.unique(df.site)


    width = 14
    height = 6
    fig = plt.figure(figsize=(width, height))
    fig.subplots_adjust(hspace=0.13)
    fig.subplots_adjust(wspace=0.13)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 16
    plt.rcParams['font.size'] = 16
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16

    labels = label_generator('lower', start="(", end=")")

    # Wombat only has one slope so don't plot it
    sites = sites[:-1]

    count=0
    for site in sites:
        df_site = df[df.site == site]
        site_name = re.sub(r"(\w)([A-Z])", r"\1 \2", site)
        ax = fig.add_subplot(2,3,1+count)

        if count == 0:
            sns.distplot(df_site.slope, ax=ax, rug=True, norm_hist=True,
                         kde_kws={"label": "KDE"})

        else:
            sns.distplot(df_site.slope, ax=ax, rug=True, norm_hist=True)

        if count < 3:
            plt.setp(ax.get_xticklabels(), visible=False)

        if count != 0 and count != 3:
            plt.setp(ax.get_yticklabels(), visible=False)

        props = dict(boxstyle='round', facecolor='white', alpha=1.0,
                     ec="white")


        fig_label = "%s %s" % (next(labels), site_name)
        ax.text(0.02, 0.95, fig_label,
                transform=ax.transAxes, fontsize=14, verticalalignment='top',
                bbox=props)
        ax.set_ylim(0, 10)
        ax.set_xlim(-0.4, 0.4)
        ax.axvline(x=0.0, ls="--", color="black", alpha=0.6)

        if count == 0:
            ax.set_ylabel("Probability density", position=(0.5, 0.0))
        if count == 4:
            #ax.set_xlabel('Temperature ($^\circ$C)', position=(1.0, 0.5))
            ax.set_xlabel('Slope (g C m$^{2}$ d$^{-1}$ $\degree$C$^{-1}$)')
        else:
            ax.set_xlabel(" ")
        from matplotlib.ticker import MaxNLocator
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.tick_params(direction='in', length=4)

        if site == "Calperum":
            ax.legend(numpoints=1, ncol=1, loc="best", frameon=False)

        count+=1


    ofdir = "/Users/mdekauwe/Dropbox/fluxnet_heatwaves_paper/figures/figs"
    ofname = "ozflux_histogram_GPP_ozflux.pdf"
    fig.savefig(os.path.join(ofdir, ofname),
                bbox_inches='tight', pad_inches=0.1)
    #plt.show()

def label_generator(case='lower', start='', end=''):
    choose_type = {'lower': string.ascii_lowercase,
                   'upper': string.ascii_uppercase}
    generator = ('%s%s%s' %(start, letter, end) for letter in choose_type[case])

    return generator

if __name__ == "__main__":

    data_dir = "data/"
    fname = os.path.join(data_dir, "ozflux_slopes_GPP.csv")
    main(fname)
