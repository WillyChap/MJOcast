import math
from pathlib import Path
import os
import pytest
import xarray as xr

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

def test_default_creation():
    '''Returns a MJO_obs instance'''
    assert create_MJOobs()

def test_default_creation2():
    '''Returns a MJO_obs_makefile instance'''
    assert create_MJOobs_file()

def test_check_obs_eofs():
    yaml_file_path = './tests/settings.yaml'
    MJO_obs = ProObs.MJOobsProcessor(yaml_file_path)
    OBS_DS, eof_list, pcs, MJO_fobs, eof_dict = MJO_obs.make_observed_MJO()
    all_pass = check_obs_eofs(MJO_obs)
    assert all_pass==True

