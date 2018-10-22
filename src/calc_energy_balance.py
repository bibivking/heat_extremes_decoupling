#!/usr/bin/env python

"""
Calculate energy balance for Ozflux sites

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (22.10.2018)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import glob
import numpy as np
import xarray as xr
import pandas as pd

import constants as c

def main(flux_dir, summary_ebr, summary_nbr, oz_flux=True):


    # wombat is not in this analysis.
    sites_we_care_about = ["Calperum","CumberlandPlains","Gingin",\
                           "GreatWesternWoodlands","Whroo","AU-ASM"]
    if oz_flux:
        flux_files = sorted(glob.glob(os.path.join(flux_dir, "*_flux.nc")))
        met_files = sorted(glob.glob(os.path.join(flux_dir, "*_met.nc")))
    else:
        #flux_files = sorted(glob.glob(os.path.join(flux_dir, "*_Flux.nc")))
        #met_files = sorted(glob.glob(os.path.join(flux_dir, "*_Met.nc")))

        # Just get Alice
        flux_files = ['data/fluxnet2015_trees/AU-ASM_2010-2013_FLUXNET2015_Flux.nc']
        met_files = ['data/fluxnet2015_trees/AU-ASM_2010-2013_FLUXNET2015_Met.nc']


    if oz_flux:
        d = get_ozflux_pfts()

    cols = ['site','pft','et','gpp_sqrt_d','gpp', 'vpd']
    df = pd.DataFrame(columns=cols)
    df2 = pd.DataFrame(columns=cols)
    for flux_fn, met_fn in zip(flux_files, met_files):

        (site, df_flx,
         df_met, pftx,
         site_name) = open_file(flux_fn, met_fn, oz_flux=oz_flux)

        if oz_flux:
            pft = d[site]
        else:
            pft = pftx

        if oz_flux:
            if pft == "EBF" or pft == "SAV":
                #print(site, np.unique(df_met.index.year))
                # daylight hours
                df_flx = df_flx.between_time("06:00", "20:00")
                df_met = df_met.between_time("06:00", "20:00")

                (df_flx, df_met) = mask_crap_days(df_flx, df_met,
                                                  oz_flux=oz_flux)
                df_met.Tair -= c.DEG_2_KELVIN

                (ebr_hw, ebr_nhw) = get_all_events(df_flx, df_met)

                if site in sites_we_care_about:
                    print(site, ebr_hw, ebr_nhw)
                    summary_ebr = np.append(summary_ebr, ebr_hw)
                    summary_nbr = np.append(summary_nbr, ebr_nhw)

        else:
            if pft == "EBF" or pft == "ENF" or pft == "DBF":

                # daylight hours
                df_flx = df_flx.between_time("06:00", "20:00")
                df_met = df_met.between_time("06:00", "20:00")

                (df_flx, df_met) = mask_crap_days(df_flx, df_met,
                                                  oz_flux=oz_flux)
                df_met.Tair -= c.DEG_2_KELVIN

                (ebr_hw, ebr_nhw) = get_all_events(df_flx, df_met)

                if site in sites_we_care_about:
                    print(site, ebr_hw, ebr_nhw)
                    summary_ebr = np.append(summary_ebr, ebr_hw)
                    summary_nbr = np.append(summary_nbr, ebr_nhw)


    return (summary_ebr, summary_nbr)

def get_all_events(df_flx, df_met):

    # Convert units ... W m-2 to kg m-2 s-1
    df_flx["ET"] = df_flx['Qle'] / latent_heat_vapourisation(df_met.Tair)

    # We need to figure out if it rained during our hot extreme as this
    # would change the Qle in the way we're searching for!
    diff = df_flx.index.minute[1] - df_flx.index.minute[0]
    # Change GPP units
    if diff == 0:
        # hour gap i.e. Tumba
        df_flx["GPP"] *= c.MOL_C_TO_GRAMS_C * c.UMOL_TO_MOL * \
                         c.SEC_TO_HR
        df_flx["ET"] *= c.SEC_TO_HR
    else:
        # 30 min gap
        df_flx["GPP"] *= c.MOL_C_TO_GRAMS_C * c.UMOL_TO_MOL * \
                         c.SEC_TO_HLFHR
        df_flx["ET"] *= c.SEC_TO_HLFHR
    tair = df_met.Tair.values
    qair = df_met.Qair.values
    press = df_met.PSurf.values
    press *= c.PA_2_KPA
    vpd = qair_to_vpd(qair, tair, press)
    #vpd = df_met.VPD.values * c.HPA_TO_KPA
    #vpd = np.where(vpd < 0.05, 0.05, vpd)

    df_flx["wue"] = df_flx["GPP"] * np.sqrt(vpd) / df_flx["ET"]
    df_dm = df_met.resample("D").max()
    df_dmu = df_met.resample("D").mean()
    df_ds = df_met.resample("D").sum()
    df_df = df_flx.resample("D").mean()
    df_dfs = df_flx.resample("D").sum()

    # We need to figure out if it rained during our hot extreme as this
    # would change the Qle in the way we're searching for!
    if diff == 0:
        # hour gap i.e. Tumba
        rain = df_met.Rainf * 3600.0
    else:
        # 30 min gap
        rain = df_met.Rainf * 1800.0
    rain = rain.fillna(0.0)
    rain = rain.resample("D").sum()

    (idx, counts, non_idx) = find_heatwave(df_dm.Tair.values, rain.values, 35.)

    if len(counts) > 0:


        df_df = df_df[~np.isnan(df_df)]

        top = df_df.Qle.values[idx] + df_df.Qh.values[idx]
        top = top[~np.isnan(top)].sum()
        #bot = np.sum(df_df.Rnet[idx] - df_df.Qg[idx])
        bot = df_df.Rnet[idx] - df_df.Qg[idx]
        bot = bot[~np.isnan(bot)].sum()
        ebr_hw =  np.mean(top / bot)

        top = df_df.Qle.values[non_idx] + df_df.Qh.values[non_idx]
        top = top[~np.isnan(top)].sum()
        #bot = np.sum(df_df.Rnet[non_idx] - df_df.Qg[non_idx])
        bot = df_df.Rnet[non_idx] - df_df.Qg[non_idx]
        #print(len(bot[~np.isnan(bot)]))
        bot = bot[~np.isnan(bot)].sum()
        ebr_nhw = np.mean(top / bot)


    else:
        ebr_hw = -999.9
        ebr_nhw = -999.9


    return (ebr_hw, ebr_nhw)


def find_heatwave(tair, rain, threshold):

    rain_threshold = 0.5 # mm d-1

    idx = np.zeros(0).astype(np.int)
    non_idx = np.zeros(0).astype(np.int)
    counts = np.zeros(0).astype(np.int)
    count = 0
    i = 0
    while i < len(tair):

        # find at least 3 consecutive days > threshold
        offset = 1

        try:
            if (tair[i] > threshold and
                tair[i+1] > threshold and
                tair[i+2] > threshold):

                while tair[i+offset] > threshold:
                    offset += 1

                end = i + offset

                # Ensure there was no rain during the heatwave and no rain 2
                # days before to exclude soil evaporation
                if (np.sum(rain[i:end]) < rain_threshold or
                    np.sum(rain[i-2:i]) < rain_threshold):
                    idx = np.append(idx, np.arange(i, end) )
                    counts = np.append(counts, np.ones(end-i) * count)

                    # count the number of heatwaves
                    count += 1
                    # we need to shift our looping index on so we skip over this
                    # heatwave
                    # -1 as the increment below will move it on one too many
                    i = end - 1
            else:
                # Exclude non-heatwave day if there was rain two days before
                # to remove soil evap contribution
                if np.sum(rain[i-2:i]) < rain_threshold:
                    non_idx = np.append(non_idx, i)

                #non_idx = np.append(non_idx, i)

        except IndexError:
            pass # trying to index too far, i.e. the next year or previous

        i += 1

    return idx, counts, non_idx

def get_ozflux_pfts():

    sites = ["AdelaideRiver","Calperum","CapeTribulation","CowBay",\
             "CumberlandPlains","DalyPasture","DalyUncleared",\
             "DryRiver","Emerald","Gingin","GreatWesternWoodlands",\
             "HowardSprings","Otway","RedDirtMelonFarm","RiggsCreek",\
             "Samford","SturtPlains","Tumbarumba","Whroo",\
             "WombatStateForest","Yanco"]

    pfts = ["SAV","EBF","TRF","TRF","EBF","GRA","SAV",\
            "SAV","NA","EBF","EBF",\
            "SAV","GRA","NA","GRA",\
            "GRA","GRA","EBF","EBF",\
            "EBF","GRA"]

    d = dict(zip(sites, pfts))

    return d

def mask_crap_days(df_flx, df_met, oz_flux=True):
    """ Mask bad QA, i.e. drop any data where Qle, Qa, Tair and Rain are flagged
    as being of poor quality"""

    if oz_flux:
        df_met.where(df_flx.Qle_qc == 1, inplace=True)
        df_met.where(df_flx.Qh_qc == 1, inplace=True)
        df_met.where(df_flx.Rnet_qc == 1, inplace=True)
        df_flx.where(df_flx.Qle_qc == 1, inplace=True)
        df_flx.where(df_flx.Qh_qc == 1, inplace=True)
        df_flx.where(df_flx.Rnet_qc == 1, inplace=True)
        df_flx.where(df_met.Tair_qc == 1, inplace=True)
        df_met.where(df_met.Tair_qc == 1, inplace=True)

    else:
        df_met.where(np.logical_or(df_flx.Qle_qc == 0, df_flx.Qle_qc == 1),
                     inplace=True)
        df_met.where(np.logical_or(df_flx.Qh_qc == 0, df_flx.Qh_qc == 1),
                     inplace=True)
        df_met.where(np.logical_or(df_flx.Rnet_qc == 0, df_flx.Rnet_qc == 1),
                     inplace=True)
        df_flx.where(np.logical_or(df_flx.Qle_qc == 0, df_flx.Qle_qc == 1),
                     inplace=True)
        df_flx.where(np.logical_or(df_flx.Qh_qc == 0, df_flx.Qh_qc == 1),
                     inplace=True)
        df_flx.where(np.logical_or(df_met.Tair_qc == 0, df_met.Tair_qc == 1),
                     inplace=True)
        df_flx.where(np.logical_or(df_flx.Rnet_qc == 0, df_flx.Rnet_qc == 1),
                     inplace=True)
        df_met.where(np.logical_or(df_met.Tair_qc == 0, df_met.Tair_qc == 1),
                     inplace=True)


    # Mask dew
    df_met.where(df_flx.Qle > 0., inplace=True)
    df_flx.where(df_flx.Qle > 0., inplace=True)

    df_met = df_met.reset_index()
    df_met = df_met.set_index('time')
    df_flx = df_flx.reset_index()
    df_flx = df_flx.set_index('time')

    return df_flx, df_met

def open_file(flux_fn, met_fn, oz_flux=True):

    if oz_flux:
        site = os.path.basename(flux_fn).split("OzFlux")[0]
    else:
        site = os.path.basename(flux_fn).split("_")[0]

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

    ds = xr.open_dataset(met_fn)
    df_met = ds.squeeze(dim=["x","y"], drop=True).to_dataframe()
    df_met = df_met.reset_index()
    df_met = df_met.set_index('time')

    return (site, df_flx, df_met, pft, site_name)

def qair_to_vpd(qair, tair, press):

    # convert back to Pa
    press /= c.PA_2_KPA

    # saturation vapor pressure
    es = 100.0 * 6.112 * np.exp((17.67 * tair) / (243.5 + tair))

    # vapor pressure
    ea = (qair * press) / (0.622 + (1.0 - 0.622) * qair)

    vpd = (es - ea) * c.PA_2_KPA
    vpd = np.where(vpd < 0.05, 0.05, vpd)

    return vpd

def latent_heat_vapourisation(tair):
    """
    Latent heat of vapourisation is approximated by a linear func of air
    temp.

    Reference:
    ----------
    * Stull, B., 1988: An Introduction to Boundary Layer Meteorology
      Boundary Conditions, pg 279.
    """
    return (2.501 - 0.00237 * tair) * 1E06

if __name__ == "__main__":

    summary_ebr = np.zeros(0)
    summary_nbr = np.zeros(0)

    print("Oz")
    flux_dir = "/Users/mdekauwe/research/OzFlux"
    (summary_ebr, summary_nbr) = main(flux_dir, summary_ebr, summary_nbr, oz_flux=True)

    # Need to get Alice Springs
    print("FLUXNET2015")
    flux_dir = "data/fluxnet2015_trees/"
    (summary_ebr, summary_nbr) = main(flux_dir, summary_ebr, summary_nbr, oz_flux=False)

    print("\n")
    print(np.min(summary_ebr), np.max(summary_ebr), np.mean(summary_ebr))
    print(np.min(summary_nbr), np.max(summary_nbr), np.mean(summary_nbr))
