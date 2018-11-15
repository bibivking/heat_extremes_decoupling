#!/usr/bin/env python

"""
For each of the OzFlux sites, plot the reln between GPP*D^0.5 and E

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

from pygam import LinearGAM
from pygam.utils import generate_X_grid

def main(fname1, fname2, fname3, fname4):

    df1 = pd.read_csv(fname1)
    df1 = df1[df1.pft == "EBF"]

    df2 = pd.read_csv(fname2)
    df2 = df2[df2.pft == "EBF"]

    # Add in Alice information from FLUXNET
    df3 = pd.read_csv(fname3)
    df3 = df3[df3.pft == "ENF"] # Alice misclassified by fluxnet
    df3 = df3[~np.isnan(df3.gpp_sqrt_d)]
    df1 = df1.append(df3, ignore_index=True)

    df4 = pd.read_csv(fname4)
    df4 = df4[df4.pft == "ENF"] # Alice misclassified by fluxnet
    df4 = df4[~np.isnan(df4.gpp_sqrt_d)]
    df2 = df2.append(df4, ignore_index=True)

    df1.loc[df1.site == "AU-ASM", "site"] = "Alice Springs"
    df2.loc[df2.site == "AU-ASM", "site"] = "Alice Springs"

    width = 14
    height = 9
    fig = plt.figure(figsize=(width, height))
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.02)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 16
    plt.rcParams['font.size'] = 16
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16

    labels = label_generator('lower', start="(", end=")")

    count=0
    sites = np.unique(df1.site)

    for site in sites:
        df_site1 = df1[df1.site == site]
        df_site2 = df2[df2.site == site]

        site_name = re.sub(r"(\w)([A-Z])", r"\1 \2", site)
        ax = fig.add_subplot(3,2,1+count)

        ax.plot(df_site2.gpp_sqrt_d, df_site2.et,
                label="Heatwave days", color="red", ls="", marker=".",
                lw=2, zorder=100, ms=9, mec="#ffb3b3")
        ax.plot(df_site1.gpp_sqrt_d, df_site1.et,
                label="Non-heatwave days", color="blue", ls="", marker=".",
                lw=2, alpha=0.3, ms=9)

        from scipy import stats

        x = df_site2.gpp_sqrt_d
        y = df_site2.et
        y = y[~np.isnan(x)]
        x = x[~np.isnan(x)]
        x = x[~np.isnan(y)]
        y = y[~np.isnan(y)]

        if site != "Alice Springs":
            nsplines = 20
        else:
            nsplines = 10

        if site != "CumberlandPlains":
            gam = LinearGAM(n_splines=nsplines).gridsearch(x, y)

            x_pred = np.linspace(min(x), max(x), num=100)
            y_pred = gam.predict(x_pred)
            y_int = gam.confidence_intervals(x_pred, width=.95)

            ax.plot(x_pred, y_pred, color="red", ls='-', lw=2.0, zorder=10)
            ax.fill_between(x_pred, y_int[:, 0], y_int[:, 1], alpha=0.2,
                            facecolor='red', zorder=10)

        x = df_site1.gpp_sqrt_d
        y = df_site1.et

        y = y[~np.isnan(x)]
        x = x[~np.isnan(x)]
        x = x[~np.isnan(y)]
        y = y[~np.isnan(y)]

        gam = LinearGAM(n_splines=nsplines).gridsearch(x, y)

        x_pred = np.linspace(min(x), max(x), num=100)
        y_pred = gam.predict(x_pred)
        y_int = gam.confidence_intervals(x_pred, width=.95)

        ax.plot(x_pred, y_pred, color="blue", ls='-', lw=2.0, zorder=10)
        ax.fill_between(x_pred, y_int[:, 0], y_int[:, 1], alpha=0.2,
                        facecolor='blue', zorder=10)

        if count < 4:
            plt.setp(ax.get_xticklabels(), visible=False)

        if count != 0 and count != 2 and count != 4:
            plt.setp(ax.get_yticklabels(), visible=False)

        props = dict(boxstyle='round', facecolor='white', alpha=1.0,
                     ec="white")


        fig_label = "%s %s" % (next(labels), site_name)
        ax.text(0.02, 0.95, fig_label,
                transform=ax.transAxes, fontsize=14, verticalalignment='top',
                bbox=props)
        ax.set_xlim(0, 11)
        ax.set_ylim(0, 4)

        if count == 2:
            ax.set_ylabel('E (mm d$^{-1}$)')
        if count == 4:
            #ax.set_xlabel('Temperature ($^\circ$C)', position=(1.0, 0.5))
            ax.set_xlabel(r"GPP $\times$ D$^{0.5}$ (g C kPa m$^{-2}$ d$^{-1}$)",
                          position=(1.0, 0.5))
        else:
            ax.set_xlabel(" ")
        from matplotlib.ticker import MaxNLocator
        ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.tick_params(direction='in', length=4)

        if site == "Calperum":
            #ax.legend(numpoints=1, ncol=1, frameon=False, loc=(1.6, 0.0))
            ax.legend(numpoints=1, ncol=1, frameon=False, loc="best")

        count+=1


    ofdir = "/Users/mdekauwe/Dropbox/fluxnet_heatwaves_paper/figures/figs"
    ofname = "ozflux_heatwave_vs_non_heatwave_gppsqrtd_et.pdf"
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
    fname1 = os.path.join(data_dir, "ozflux_heatwaves_gpp_sqrt_d_et_non.csv")
    fname2 = os.path.join(data_dir, "ozflux_heatwaves_gpp_sqrt_d_et.csv")
    fname3 = os.path.join(data_dir,
                         "fluxnet2015_heatwaves_gpp_sqrt_d_et_non.csv")
    fname4 = os.path.join(data_dir, "fluxnet2015_heatwaves_gpp_sqrt_d_et.csv")
    main(fname1, fname2, fname3, fname4)
