import math
from pathlib import Path
import os
import pytest
import xarray as xr
import glob

import numpy as np
#import pytest
import MJOcast.utils.WHtools as whtools 
import MJOcast.utils.ProcessForecasts as ProFo 
import MJOcast.utils.ProcessOBS as ProObs
# FIXME: Define paths as fixture at a central point

def create_MJOobs():
    print(os.getcwd())
    yaml_file_path = './tests/settings.yaml'
    MJO_obs = ProObs.MJOobsProcessor(yaml_file_path)
    return MJO_obs

def create_MJOobs_file():
    yaml_file_path = './tests/settings.yaml'
    MJO_obs = ProObs.MJOobsProcessor(yaml_file_path)
    OBS_DS, eof_list, pcs, MJO_fobs, eof_dict = MJO_obs.make_observed_MJO()
    return MJO_fobs

def create_For_file_case1():
    yaml_file_path = './tests/settings.yaml'
    MJO_obs = ProObs.MJOobsProcessor(yaml_file_path)
    OBS_DS, eof_list, pcs, MJO_fobs, eof_dict = MJO_obs.make_observed_MJO()
    MJO_for = ProFo.MJOforecaster(yaml_file_path,MJO_obs.eof_dict,MJO_obs.MJO_fobs)
    
    svname = glob.glob('./tests//MJO_Forecast_Init_Case*.nc')
    for sv in svname: 
        if os.path.exists(sv):
            os.remove(sv)
        
    DS_CESM_for,OLR_cesm_anom_filterd,U200_cesm_anom_filterd,U850_cesm_anom_filterd = MJO_for.create_forecasts()
    
    all_pass = corr_rmm1_rmm2(MJO_for)
    
    svname = glob.glob('./tests//MJO_Forecast_Init_Case*.nc')
    for sv in svname: 
        if os.path.exists(sv):
            os.remove(sv)
    
    return all_pass

def create_For_file_case2():
    yaml_file_path = './tests/settings_LatLonShift_newlons.yaml'
    MJO_obs = ProObs.MJOobsProcessor(yaml_file_path)
    OBS_DS, eof_list, pcs, MJO_fobs, eof_dict = MJO_obs.make_observed_MJO()
    MJO_for = ProFo.MJOforecaster(yaml_file_path,MJO_obs.eof_dict,MJO_obs.MJO_fobs)
    
    svname = glob.glob('./tests//MJO_Forecast_Init_Case*.nc')
    for sv in svname: 
        if os.path.exists(sv):
            os.remove(sv)
        
    
    DS_CESM_for,OLR_cesm_anom_filterd,U200_cesm_anom_filterd,U850_cesm_anom_filterd = MJO_for.create_forecasts()
    
    all_pass = corr_rmm1_rmm2(MJO_for)
    
    svname = glob.glob('./tests//MJO_Forecast_Init_Case*.nc')
    for sv in svname: 
        if os.path.exists(sv):
            os.remove(sv)
    
    
    return all_pass

def corr_rmm1_rmm2(MJO_for):
    """
    Verify forecast RMM1/RMM2 match the reference in shape AND magnitude.

    A correlation > 0.7 alone is scale-invariant and would miss a scaling error, so we
    tighten the correlation (catches sign/structure) and check the standard-deviation
    ratio (catches scaling). Absolute values are not pinned because the EOF solve differs
    slightly across environments (~0.3 in RMM units).
    """
    fp = './tests/test_cases/MJO_Forecast_Init_standard.nc'
    check_corr_ds = xr.open_dataset(fp)

    check_modes = ['RMM1', 'RMM2']

    for cm in check_modes:
        ref = check_corr_ds[cm].sel(number=1).isel(time=slice(0,12)).values.squeeze()
        if len(MJO_for.MJO_forecast_DS.number) > 1:
            got = MJO_for.MJO_forecast_DS[cm].sel(number=1).isel(time=slice(0,12)).values.squeeze()
        else:
            got = MJO_for.MJO_forecast_DS[cm].isel(time=slice(0,12)).values.squeeze()
        corr = np.corrcoef(ref, got)[0, 1]
        std_ratio = np.std(got) / np.std(ref)
        assert corr > 0.9, f"{cm}: correlation {corr:.3f} <= 0.9 (sign/structure error)"
        assert 0.75 < std_ratio < 1.3, f"{cm}: std ratio {std_ratio:.3f} (scaling error)"

    return True

def eof1_u200_matches_obs():
    """Regression check: the saved forecast eof1_u200 must be EOF1 (not a copy of eof2_u200)."""
    yaml_file_path = './tests/settings.yaml'
    MJO_obs = ProObs.MJOobsProcessor(yaml_file_path)
    MJO_obs.make_observed_MJO()
    MJO_for = ProFo.MJOforecaster(yaml_file_path, MJO_obs.eof_dict, MJO_obs.MJO_fobs)

    for sv in glob.glob('./tests//MJO_Forecast_Init_Case*.nc'):
        if os.path.exists(sv):
            os.remove(sv)

    MJO_for.create_forecasts()
    ds = MJO_for.MJO_forecast_DS

    distinct = not np.array_equal(ds['eof1_u200'].values, ds['eof2_u200'].values)
    matches_obs = np.array_equal(ds['eof1_u200'].values, np.array(MJO_obs.MJO_fobs['eof1_u200']))

    for sv in glob.glob('./tests//MJO_Forecast_Init_Case*.nc'):
        if os.path.exists(sv):
            os.remove(sv)

    return distinct and matches_obs


def test_eof1_u200_matches_obs():
    '''eof1_u200 in the forecast output must be EOF1, not a duplicate of eof2_u200'''
    assert eof1_u200_matches_obs() == True


def test_default_creation_case1():
    '''Returns a MJO_obs instance'''
    all_pass = create_For_file_case1()
    assert all_pass==True
    

def test_default_creation_case2():
    '''Returns a MJO_obs instance'''
    all_pass = create_For_file_case2()
    assert all_pass==True
    


def check_obs_eofs(MJO_obs):
    """
    Checks the correlation between observed dataset EOF values and ERA5 EOF values for specific modes.

    This function opens an observed dataset file, interpolates the data, and calculates the correlation coefficients
    between specific EOF modes of the observed dataset and corresponding ERA5 EOF values. It compares the calculated
    correlations to a predefined threshold and raises an assertion error if any correlation is below the threshold.

    Returns:
        bool: True if all correlations are above the threshold, indicating a high correlation between observed and ERA5 EOF values.
    """
    fp = './tests/test_cases/eofs_MJO.nc'
    check_corr_ds = xr.open_dataset(fp)
    check_corr_ds = whtools.interpolate_obs(check_corr_ds, MJO_obs.forecast_lons)
    corr_dict = {}

    check_modes = ['eof1_olr', 'eof2_olr', 'eof1_u200', 'eof2_u200', 'eof1_u850', 'eof2_u850']

    for cm in check_modes:
        corrnum = np.corrcoef(check_corr_ds[cm].values, MJO_obs.MJO_fobs[cm])[0, 1]
        corr_dict[cm] = corrnum

    print('the correlation of the observed dataset EOF values and the ERA5 EOF values is:')
    print(corr_dict)

    for cm in check_modes:
        
        if corr_dict[cm] > 0.7:
            all_pass=True
        else:
            all_pass=False
            break

    return all_pass ==True

