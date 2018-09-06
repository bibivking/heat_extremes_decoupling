#!/usr/bin/env python

"""
Get the years for the tabe in the paper...

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (01.09.2018)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import numpy as np
import xarray as xr
import pandas as pd

import constants as c

def main(flux_dir, ofname, oz_flux=True):

    if oz_flux:
        sites = ["Calperum","CumberlandPlains","Gingin",\
                 "GreatWesternWoodlands","Whroo","WombatStateForest","Yanco"]
    else:
        sites = ["AU-ASM","IT-CA1","IT-CA3","FR-LBr","FR-Pue",\
                 "ZM-Mon","IT-PT1","CN-Qia","IT-Ro1", "IT-Ro2", "US-MMS"]

    for site in sites:
        if oz_flux:
            flux_fn = os.path.join(flux_dir, "%sOzFlux2.0_flux.nc" % (site))
        else:
            flux_fn = glob.glob(os.path.join(flux_dir,
                                "%s_*_FLUXNET2015_Flux.nc" % (site)))[0]

        ds = xr.open_dataset(flux_fn)

        if oz_flux == False:
            pft = (str(ds.IGBP_veg_short.values, 'utf-8'))
            pft = pft.replace(" ", "")

            site_name = ds.site_name
            site_name = site_name.replace(" ", "")
        else:
            pft = None
            site_name = None
        df_flx = ds.squeeze(dim=["x","y"], drop=True).to_dataframe()
        df_flx = df_flx.reset_index()
        df_flx = df_flx.set_index('time')

        print(site, np.unique(df_flx.index.year))


if __name__ == "__main__":

    print("Oz")
    flux_dir = "/Users/mdekauwe/research/OzFlux"
    ofname = "ozflux_all_events.csv"
    main(flux_dir, ofname, oz_flux=True)
    print("\n")
    print("FLUXNET2015")
    flux_dir = "data/fluxnet2015_trees"
    ofname = "fluxnet2015_all_events.csv"
    main(flux_dir, ofname, oz_flux=False)
