#!/usr/bin/env python

"""
For each of the FLUXNET2015 sites, plot the TXx and T-3 days
GPP slope as a function of mean daily air temp

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (20.04.2018)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import netCDF4 as nc
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import re
import constants as c
import string

def main(fname, slope_type):

    df = pd.read_csv(fname)
    df = df[~np.isnan(df.temp)]

    # Drop disturbance sites
    ignore_sites = ["UCI-1850burnsite", "UCI-1964burnsitewet", \
                    "Saskatchewan-WesternBoreal,forestburnedin1989"]
    for site in ignore_sites:
        df = df.drop( df[(df.site == site)].index )

    # Drop repeated Aus sites and Loxton which is an irrigated orchard
    ignore_sites = ["AliceSprings",  "WallabyCreek",\
                    "Whroo", "CumberlandPlains", "Loxton"]
    for site in ignore_sites:
        df = df.drop( df[(df.site == site)].index )

    # no data for matching gpp and GPP slopes, so drop from dataframe so we
    # don't plot empty panels
    ignore_sites = ["ParcoTicinoforest", "Virasoro"]
    for site in ignore_sites:
        df = df.drop( df[(df.site == site)].index )

    # clean names
    df.loc[df.site == "Casteld'Asso1", "site"] = "Castel d'Asso"
    df.loc[df.site == "Casteld'Asso3", "site"] = "Castel d'Asso"
    df.loc[df.site == "Roccarespampani1", "site"] = "Roccarespampani"
    df.loc[df.site == "Roccarespampani2", "site"] = "Roccarespampani"
    df.loc[df.site == "ParcoTicinoforest", "site"] = "Parco Ticino forest"

    #width  = 12.0
    #height = width / 1.618
    #print(width, height)
    #sys.exit()
    width = 10
    height = 8
    fig = plt.figure(figsize=(width, height))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    labels = label_generator('lower', start="(", end=")")

    count = 0
    sites = np.unique(df.site)
    for site in sites:
        site_name = re.sub(r"(\w)([A-Z])", r"\1 \2", site)
        print(site)

        if site == "Casteld'Asso1":
            site_name = "Castel d'Asso1"
        elif site == "Casteld'Asso3":
            site_name = "Castel d'Asso3"
        elif site == "Casteld'Asso":
            site_name = "Castel d'Asso"
        #site_name = site[:6]

        ax = fig.add_subplot(4,2,1+count)

        df_site = df[df.site == site]
        events = int(len(df_site)/4)

        cnt = 0
        for e in range(0, events):

            from scipy import stats
            x = df_site["temp"][cnt:cnt+4]
            y = df_site["GPP"][cnt:cnt+4]
            (slope, intercept,
             r_value, p_value, std_err) = stats.linregress(x,y)

            # We only want to plot slopes if the direction of change in GPP is
            # opposite the LE, i.e. that they change together
            xx = df_site["temp"][cnt:cnt+4]
            yy = df_site["Qle"][cnt:cnt+4]
            (slope2, intercept2, r_value2,
             p_value2, std_err2) = stats.linregress(xx,yy)

            if slope_type == "negative":
                if (slope < 0.0 and p_value <= 0.05) and (slope2 > 0.0):
                    ax.plot(df_site["temp"][cnt:cnt+4],
                            df_site["GPP"][cnt:cnt+4],
                            label=site, color="darkblue", ls="-", marker="o",
                            zorder=100, alpha=0.5, lw=2.0)
                elif (slope < 0.0 and p_value > 0.05) and (slope2 > 0.0):
                    ax.plot(df_site["temp"][cnt:cnt+4],
                            df_site["GPP"][cnt:cnt+4],
                            label=site, ls="-", marker="o", color="darkgreen",
                            alpha=0.4, zorder=1)
            elif slope_type == "positive":
                if (slope > 0.0 and p_value <= 0.05) and (slope2 > 0.0):
                    ax.plot(df_site["temp"][cnt:cnt+4],
                            df_site["GPP"][cnt:cnt+4],
                            label=site, color="darkblue", ls="-", marker="o",
                            zorder=100, alpha=0.5, lw=2.0)
                elif (slope > 0.0 and p_value > 0.05) and (slope2 > 0.0):
                    ax.plot(df_site["temp"][cnt:cnt+4],
                            df_site["GPP"][cnt:cnt+4],
                            label=site, ls="-", marker="o", color="darkgreen",
                            zorder=1, alpha=0.4)

            cnt += 4

        if count == 4:
            ax.set_ylabel("GPP (g C m$^{-2}$ d$^{-1}$)", position=(0.5, 1.0))
        if count == 6:
            ax.set_xlabel('Temperature ($\degree$C)', position=(1.0, 0.5))
            #ax.set_xlabel('Temperature ($^\circ$C)')

        if count < 5:
            plt.setp(ax.get_xticklabels(), visible=False)

        if count != 0 and count != 2 and count != 4 and count != 6:
            plt.setp(ax.get_yticklabels(), visible=False)

        props = dict(boxstyle='round', facecolor='white', alpha=1.0,
                     ec="white")
        fig_label = "%s %s" % (next(labels), site_name)
        ax.text(0.02, 0.95, fig_label,
                transform=ax.transAxes, fontsize=12, verticalalignment='top',
                bbox=props)

        from matplotlib.ticker import MaxNLocator
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.set_ylim(0, 15)
        ax.set_xlim(18, 47)
        count += 1

    ofdir = "/Users/mdekauwe/Dropbox/fluxnet_heatwaves_paper/figures/figs"
    ofname = "all_events_GPP_FLUXNET_%s.pdf" % (slope_type)
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
    fname = "fluxnet2015_all_events.csv"
    fname = os.path.join(data_dir, fname)
    main(fname, slope_type="negative")
    main(fname, slope_type="positive")
