#!/usr/bin/env python

"""
For each extreme event (i.e. where the hottest day > 37 degC and the three days
prior), record the Qle, GPP, Bowen ratio.

- We are screening "bad" flux data
- We are dropping events that have gaps
- We are dropping events where it rained as this would mess with the Qle.

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (26.04.2017)"
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
        flux_files = sorted(glob.glob(os.path.join(flux_dir, "*_flux.nc")))
        met_files = sorted(glob.glob(os.path.join(flux_dir, "*_met.nc")))
    else:
        flux_files = sorted(glob.glob(os.path.join(flux_dir, "*_Flux.nc")))
        met_files = sorted(glob.glob(os.path.join(flux_dir, "*_Met.nc")))

    output_dir = "data"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if oz_flux:
        d = get_ozflux_pfts()

    cols = ['site','pft','temp','Qle','GPP','vpd','wue','yr','doy']
    df = pd.DataFrame(columns=cols)
    for flux_fn, met_fn in zip(flux_files, met_files):
        (site, df_flx,
         df_met, pftx,
         site_name) = open_file(flux_fn, met_fn, oz_flux=oz_flux)
        print(site)
        if oz_flux:
            pft = d[site]
        else:
            pft = pftx

        if oz_flux:
            # We aren't using SAV sites, but pass them through for ozflux files
            if pft == "EBF" or pft == "SAV" or pft == "TRF":
                #print(site, np.unique(df_met.index.year))

                (df_flx, df_met) = mask_crap_days(df_flx, df_met,
                                                  oz_flux=oz_flux)

                (Tairs, Qles, GPPs, vpds, wues,
                idx_yrs, idx_doy) = get_all_events(df_flx, df_met)

                lst = []
                for i in range(len(Tairs)):
                    lst.append([site,pft,Tairs[i],Qles[i],GPPs[i],\
                                vpds[i],wues[i],idx_yrs[i],idx_doy[i]])
                dfx = pd.DataFrame(lst, columns=cols)
                # reverse the order hot to cool
                #dfx = dfx.reindex(index=dfx.index[::-1])
                df = df.append(dfx)
        else:
            if pft == "EBF" or pft == "ENF" or pft == "DBF":

                (df_flx, df_met) = mask_crap_days(df_flx, df_met,
                                                  oz_flux=oz_flux)

                (Tairs, Qles, GPPs, vpds, wues,
                idx_yrs, idx_doy) = get_all_events(df_flx, df_met)

                lst = []
                for i in range(len(Tairs)):
                    lst.append([site_name,pft,Tairs[i],Qles[i],GPPs[i],\
                                vpds[i],wues[i],idx_yrs[i],idx_doy[i]])

                dfx = pd.DataFrame(lst, columns=cols)
                # reverse the order hot to cool
                #dfx = dfx.reindex(index=dfx.index[::-1])
                df = df.append(dfx)

    df.to_csv(os.path.join(output_dir, ofname), index=False)

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
    press *= c.PA_TO_KPA
    vpd = qair_to_vpd(qair, tair, press)

    df_flx["wue"] = df_flx["GPP"] * np.sqrt(vpd) / df_flx["ET"]

    df_dm = df_met.resample("D").max()
    df_dmu = df_met.resample("D").mean()
    df_ds = df_met.resample("D").sum()
    df_df = df_flx.resample("D").mean()
    df_dfs = df_flx.resample("D").sum()

    # We need to figure out if it rained during our hot extreme as this
    # would change the Qle in the way we're searching for!
    diff = df_met.index.minute[1] - df_met.index.minute[0]
    if diff == 0:
        # hour gap i.e. Tumba
        rain = df_met.Rainf * 3600.0
    else:
        # 30 min gap
        rain = df_met.Rainf * 1800.0
    rain = rain.fillna(0.0)
    rain = rain.resample("D").sum()

    #
    ## Get the TXx event
    #
    TXx = df_dm.sort_values("Tair", ascending=False)[:1].Tair.values[0]
    TXx_idx = df_dm.sort_values("Tair", ascending=False)[:1].index.values[0]
    TXx_idx_minus_four = TXx_idx - pd.Timedelta(3, unit='d')

    # Need this to screen for rain in the 48 hours before the event, i.e.
    # to remove soil evap contributions
    TXx_idx_minus_six = TXx_idx - pd.Timedelta(5, unit='d')

    (Tairs, Qles,
     Qhs, GPP,
     vpds, wues,
     idx_yrs, idx_doy) = get_values(df_dm, df_dmu, df_df, df_dfs, TXx_idx,
                                    TXx_idx_minus_four, TXx_idx_minus_six)

    (Tairs, Qles, GPP, vpds, wues,
     df_dm, df_dmu, df_df,
     idx_yrs, idx_doy) = is_event_long_enough(df_dm,df_dmu, df_df, df_dfs,
                                              TXx_idx, TXx_idx_minus_four,
                                              Tairs, Qles, GPP, rain, vpds,
                                              wues, idx_yrs, idx_doy)

    if len(Tairs) < 4:
        Tairs = np.array([np.nan,np.nan,np.nan,np.nan])
        Qles = np.array([np.nan,np.nan,np.nan,np.nan])
        GPP = np.array([np.nan,np.nan,np.nan,np.nan])
        vpds = np.array([np.nan,np.nan,np.nan,np.nan])
        wues = np.array([np.nan,np.nan,np.nan,np.nan])
        idx_yrs = np.array([np.nan,np.nan,np.nan,np.nan])
        idx_doy = np.array([np.nan,np.nan,np.nan,np.nan])

    if (TXx) < 37.0:
        Tairs = np.array([np.nan,np.nan,np.nan,np.nan])
        Qles = np.array([np.nan,np.nan,np.nan,np.nan])
        GPP = np.array([np.nan,np.nan,np.nan,np.nan])
        vpds = np.array([np.nan,np.nan,np.nan,np.nan])
        wues = np.array([np.nan,np.nan,np.nan,np.nan])
        idx_yrs = np.array([np.nan,np.nan,np.nan,np.nan])
        idx_doy = np.array([np.nan,np.nan,np.nan,np.nan])

    #
    ## Get all the events other than the TXx that are > Tthresh
    #

    # Drop the hottest event as we've already got it
    df_dmu = df_dmu[(df_dm.index < TXx_idx_minus_four) |
                   (df_dm.index > TXx_idx)]
    df_dm = df_dm[(df_dm.index < TXx_idx_minus_four) |
                   (df_dm.index > TXx_idx)]
    df_df = df_df[(df_df.index < TXx_idx_minus_four) |
                   (df_df.index > TXx_idx)]

    # Then get next TXx
    TXx = df_dm.sort_values("Tair", ascending=False)[:1].Tair.values[0]

    while TXx > 37.0:

        # Then get next TXx
        TXx = df_dm.sort_values("Tair", ascending=False)[:1].Tair.values[0]
        TXx_idx = df_dm.sort_values("Tair", ascending=False)[:1].index.values[0]
        TXx_idx_minus_four = TXx_idx - pd.Timedelta(3, unit='d')

        # Need this to screen for rain in the 48 hours before the event, i.e.
        # to remove soil evap contributions
        TXx_idx_minus_six = TXx_idx - pd.Timedelta(5, unit='d')

        (Tairsx, Qlesx,
         Qhsx, GPPx,
         vpdsx, wuesx,
         idx_yrsx, idx_doyx) = get_values(df_dm, df_dmu, df_df, df_dfs, TXx_idx,
                                          TXx_idx_minus_four, TXx_idx_minus_six)

        (Tairsx, Qlesx,
         GPPx, vpdsx, wuesx,
         df_dm, df_dmu, df_df,
         idx_yrsx, idx_doyx) = is_event_long_enough(df_dm, df_dmu, df_df,
                                                    df_dfs, TXx_idx,
                                                    TXx_idx_minus_four,
                                                    Tairsx, Qlesx, GPPx,
                                                    rain, vpdsx, wuesx,
                                                    idx_yrsx, idx_doyx)

        Tairsx = Tairsx[~np.isnan(Tairsx)]
        Qlesx = Qlesx[~np.isnan(Qlesx)]
        GPPx = GPPx[~np.isnan(GPPx)]
        vpdsx = vpdsx[~np.isnan(vpdsx)]
        wuesx = wuesx[~np.isnan(wuesx)]
        idx_yrsx = idx_yrsx[~np.isnan(idx_yrsx)]
        idx_doyx = idx_doyx[~np.isnan(idx_doyx)]

        if len(Tairsx) == 4:
            Tairs = np.append(Tairs, Tairsx)
            Qles = np.append(Qles, Qlesx)
            GPP = np.append(GPP, GPPx)
            vpds = np.append(vpds, vpdsx)
            wues = np.append(wues, wuesx)
            idx_yrs = np.append(idx_yrs, idx_yrsx)
            idx_doy = np.append(idx_doy, idx_doyx)

        # Drop this event
        df_dmu = df_dmu[(df_dm.index < TXx_idx_minus_four) |
                       (df_dm.index > TXx_idx)]
        df_dm = df_dm[(df_dm.index < TXx_idx_minus_four) |
                       (df_dm.index > TXx_idx)]
        df_df = df_df[(df_df.index < TXx_idx_minus_four) |
                       (df_df.index > TXx_idx)]

        rain = rain[(rain.index < TXx_idx_minus_four) |
                    (rain.index > TXx_idx)]

        # Then get next TXx
        TXx = df_dm.sort_values("Tair",
                                ascending=False)[:1].Tair.values[0]

    Tairs = Tairs[~np.isnan(Tairs)]
    Qles = Qles[~np.isnan(Qles)]
    GPP = GPP[~np.isnan(GPP)]
    vpds = vpds[~np.isnan(vpds)]
    wues = wues[~np.isnan(wues)]
    idx_yrs = idx_yrs[~np.isnan(idx_yrs)]
    idx_doy = idx_doy[~np.isnan(idx_doy)]

    if len(Tairs) < 4:
        Tairs = np.array([np.nan,np.nan,np.nan,np.nan])
        Qles = np.array([np.nan,np.nan,np.nan,np.nan])
        GPP = np.array([np.nan,np.nan,np.nan,np.nan])
        vpds = np.array([np.nan,np.nan,np.nan,np.nan])
        wues = np.array([np.nan,np.nan,np.nan,np.nan])
        idx_yrs = np.array([np.nan,np.nan,np.nan,np.nan])
        idx_doy = np.array([np.nan,np.nan,np.nan,np.nan])

    return(Tairs, Qles, GPP, vpds, wues, idx_yrs, idx_doy)

def get_values(df_dm, df_dmu, df_df, df_dfs, TXx_idx, TXx_idx_minus_four,
               TXx_idx_minus_six):

    Tairs = df_dm[(df_dm.index >= TXx_idx_minus_four) &
                  (df_dm.index <= TXx_idx)].Tair.values
    Qles = df_df[(df_dm.index >= TXx_idx_minus_four) &
                 (df_dm.index <= TXx_idx)].Qle.values
    Qhs = df_df[(df_dm.index >= TXx_idx_minus_four) &
                (df_dm.index <= TXx_idx)].Qh.values

    GPPs = df_dfs[(df_dfs.index >= TXx_idx_minus_four) &
                  (df_dfs.index <= TXx_idx)].GPP.values


    qair = df_dmu[(df_dm.index >= TXx_idx_minus_four) &
                  (df_dm.index <= TXx_idx)].Qair.values
    press = df_dmu[(df_dm.index >= TXx_idx_minus_four) &
                  (df_dm.index <= TXx_idx)].PSurf.values
    press *= c.PA_TO_KPA


    vpds = qair_to_vpd(qair, Tairs, press)

    wues = df_df[(df_dm.index >= TXx_idx_minus_four) &
                 (df_dm.index <= TXx_idx)].wue.values


    idx_yrs = df_dm[(df_dm.index >= TXx_idx_minus_four) &
                    (df_dm.index <= TXx_idx)].index.year.values
    idx_doy = df_dm[(df_dm.index >= TXx_idx_minus_four) &
                    (df_dm.index <= TXx_idx)].index.dayofyear.values

    return (Tairs, Qles, Qhs, GPPs, vpds, wues, idx_yrs, idx_doy)

def is_event_long_enough(df_dm, df_dmu, df_df, df_dfs, TXx_idx,
                         TXx_idx_minus_four, Tairs, Qles, GPPs, rain, vpds,
                         wues, idx_yrs, idx_doy):

    while len(Tairs) != 4:

        # Drop this event as it wasn't long enough
        df_dmu = df_dmu[(df_dm.index < TXx_idx_minus_four) |
                       (df_dm.index > TXx_idx)]
        df_dm = df_dm[(df_dm.index < TXx_idx_minus_four) |
                       (df_dm.index > TXx_idx)]
        df_df = df_df[(df_df.index < TXx_idx_minus_four) |
                       (df_df.index > TXx_idx)]


        TXx = df_dm.sort_values("Tair", ascending=False)[:1].Tair.values[0]
        TXx_idx = df_dm.sort_values("Tair",
                                    ascending=False)[:1].index.values[0]
        TXx_idx_minus_four = TXx_idx - pd.Timedelta(3, unit='d')

        # Need this to screen for rain in the 48 hours before the event, i.e.
        # to remove soil evap contributions
        TXx_idx_minus_six = TXx_idx - pd.Timedelta(5, unit='d')

        (Tairs, Qles,
         Qhs, GPPs,
         vpds, wues,
         idx_yrs, idx_doy) = get_values(df_dm, df_dmu, df_df, df_dfs, TXx_idx,
                                        TXx_idx_minus_four, TXx_idx_minus_six)

        (Tairs, Qles, GPPs, vpds, wues,
         df_dm, df_dmu, df_df,
         idx_yrs, idx_doy) = check_for_rain(rain, TXx_idx_minus_four,
                                            TXx_idx_minus_six,
                                            TXx_idx, df_dm, df_dmu, df_df,
                                            df_dfs, Tairs, Qles, Qhs, GPPs,
                                            vpds, wues, idx_yrs, idx_doy)

        if len(df_dm) <= 4:
            Tairs = np.array([np.nan,np.nan,np.nan,np.nan])
            Qles = np.array([np.nan,np.nan,np.nan,np.nan])
            GPP = np.array([np.nan,np.nan,np.nan,np.nan])
            vpds = np.arrary([np.nan,np.nan,np.nan,np.nan])
            wues = np.arrary([np.nan,np.nan,np.nan,np.nan])
            idx_yrs = np.arrary([np.nan,np.nan,np.nan,np.nan])
            idx_doy = np.arrary([np.nan,np.nan,np.nan,np.nan])
            break

    return (Tairs, Qles, GPPs, vpds, wues, df_dm, df_dmu, df_df,
            idx_yrs, idx_doy)



def check_for_rain(rain, TXx_idx_minus_four, TXx_idx_minus_six,
                   TXx_idx, df_dm, df_dmu, df_df,
                   df_dfs, Tairs, Qles, Qhs, GPPs, vpds, wues,
                   idx_yrs, idx_doy):

    threshold = 0.5 # mm d-1
    total_rain = np.sum(rain[(rain.index >= TXx_idx_minus_four) &
                             (rain.index <= TXx_idx)].values)

    # rain 2 days before, i.e. rule out soil evap contributions
    total_rain_prior = np.sum(rain[(rain.index >= TXx_idx_minus_six) &
                                   (rain.index <= TXx_idx_minus_four)].values)

    while total_rain > threshold or total_rain_prior > threshold or \
           len(Tairs) != 4:
    #while total_rain > threshold  or len(Tairs) != 4:

        # Drop this event as there was some rain or we didn't get 5 good QA days
        df_dmu = df_dmu[(df_dm.index < TXx_idx_minus_four) |
                       (df_dm.index > TXx_idx)]
        df_dm = df_dm[(df_dm.index < TXx_idx_minus_four) |
                       (df_dm.index > TXx_idx)]
        df_df = df_df[(df_df.index < TXx_idx_minus_four) |
                       (df_df.index > TXx_idx)]

        # Get a new event
        TXx = df_dm.sort_values("Tair", ascending=False)[:1].Tair.values[0]
        TXx_idx = df_dm.sort_values("Tair",
                                    ascending=False)[:1].index.values[0]
        TXx_idx_minus_four = TXx_idx - pd.Timedelta(3, unit='d')

        # Need this to screen for rain in the 48 hours before the event, i.e.
        # to remove soil evap contributions
        TXx_idx_minus_six = TXx_idx - pd.Timedelta(5, unit='d')

        (Tairs, Qles,
         Qhs, GPPs,
         vpds, wues,
         idx_yrs, idx_doy) = get_values(df_dm, df_dmu, df_df, df_dfs, TXx_idx,
                                        TXx_idx_minus_four, TXx_idx_minus_six)

        total_rain = np.sum(rain[(rain.index >= TXx_idx_minus_four) &
                                 (rain.index <= TXx_idx)].values)

        # rain 2 days before, i.e. rule out soil evap contributions
        total_rain_prior = np.sum(rain[(rain.index >= TXx_idx_minus_six) &
                                  (rain.index <= TXx_idx_minus_four)].values)


        if len(df_dm) <= 4:
            Tairs = np.array([np.nan,np.nan,np.nan,np.nan])
            Qles = np.array([np.nan,np.nan,np.nan,np.nan])
            GPPs = np.array([np.nan,np.nan,np.nan,np.nan])
            vpds = np.array([np.nan,np.nan,np.nan,np.nan])
            wues = np.array([np.nan,np.nan,np.nan,np.nan])
            idx_yrs = np.array([np.nan,np.nan,np.nan,np.nan])
            idx_doy = np.array([np.nan,np.nan,np.nan,np.nan])
            break

    return (Tairs, Qles, GPPs, vpds, wues, df_dm, df_dmu, df_df, idx_yrs,
            idx_doy)



def get_ozflux_pfts():

    sites = ["AdelaideRiver","Calperum","CapeTribulation","CowBay",\
             "CumberlandPlains","DalyPasture","DalyUncleared",\
             "DryRiver","Emerald","Gingin","GreatWesternWoodlands",\
             "HowardSprings","Otway","RedDirtMelonFarm","RiggsCreek",\
             "Samford","SturtPlains","Tumbarumba","Whroo",\
             "WombatStateForest","Yanco"]

    pfts = ["SAV","EBF","TRF","TRF",\
            "EBF","GRA","SAV",\
            "SAV","NA","EBF","EBF",\
            "SAV","GRA","NA","GRA",\
            "GRA","GRA","EBF","EBF",\
            "EBF","GRA"]

    d = dict(zip(sites, pfts))

    return d

def mask_crap_days(df_flx, df_met, oz_flux=True):
    """ Mask bad QA, i.e. drop any data where Qle, Qa, Tair and Rain are flagged
    as being of poor quality"""

    # daylight hours
    df_flx = df_flx.between_time("06:00", "20:00")
    df_met = df_met.between_time("06:00", "20:00")

    if oz_flux:
        df_met.where(df_flx.Qle_qc == 1, inplace=True)
        #df_met.where(df_flx.Qh_qc == 1, inplace=True)
        df_flx.where(df_flx.Qle_qc == 1, inplace=True)
        #df_flx.where(df_flx.Qh_qc == 1, inplace=True)
        df_flx.where(df_met.Tair_qc == 1, inplace=True)
        df_met.where(df_met.Tair_qc == 1, inplace=True)

    else:
        df_met.where(np.logical_or(df_flx.Qle_qc == 0, df_flx.Qle_qc == 1),
                     inplace=True)
        #df_met.where(np.logical_or(df_flx.Qh_qc == 0, df_flx.Qh_qc == 1),
        #             inplace=True)

        df_flx.where(np.logical_or(df_flx.Qle_qc == 0, df_flx.Qle_qc == 1),
                     inplace=True)
        #df_flx.where(np.logical_or(df_flx.Qh_qc == 0, df_flx.Qh_qc == 1),
        #             inplace=True)
        df_flx.where(np.logical_or(df_met.Tair_qc == 0, df_met.Tair_qc == 1),
                     inplace=True)

        df_met.where(np.logical_or(df_met.Tair_qc == 0, df_met.Tair_qc == 1),
                     inplace=True)


    # Mask dew
    df_met.where(df_flx.Qle > 0., inplace=True)
    df_flx.where(df_flx.Qle > 0., inplace=True)

    df_met.Tair -= c.DEG_2_KELVIN

    df_met = df_met.reset_index()
    df_met = df_met.set_index('time')
    df_flx = df_flx.reset_index()
    df_flx = df_flx.set_index('time')

    return df_flx, df_met

def open_file(flux_fn, met_fn, oz_flux=True):
    site = os.path.basename(flux_fn).split("OzFlux")[0]
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
    press /= c.PA_TO_KPA

    # saturation vapor pressure
    es = 100.0 * 6.112 * np.exp((17.67 * tair) / (243.5 + tair))

    # vapor pressure
    ea = (qair * press) / (0.622 + (1.0 - 0.622) * qair)

    vpd = (es - ea) * c.PA_TO_KPA
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

    print("Oz")
    flux_dir = "/Users/mdekauwe/research/OzFlux"
    ofname = "ozflux_all_events.csv"
    main(flux_dir, ofname, oz_flux=True)

    print("FLUXNET2015")
    flux_dir = "/Users/mdekauwe/research/heat_extremes_decoupling/data/fluxnet2015_trees"
    ofname = "fluxnet2015_all_events.csv"
    main(flux_dir, ofname, oz_flux=False)
