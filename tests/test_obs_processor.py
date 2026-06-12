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
    Verify each observed EOF mode matches the ERA5 reference in shape AND magnitude.

    A bare correlation > 0.7 would pass even with a sign flip masked by structure, or a
    2x scaling error (correlation is scale-invariant). So in addition to a much tighter
    correlation we check the standard-deviation ratio (catches scaling) and the absolute
    values within a tolerance that still accommodates EOF-solver/LAPACK differences
    across environments (measured ~1e-2 on values spanning ~0.13).

    Returns:
        bool: True if every mode passes; raises AssertionError naming the failing mode otherwise.
    """
    fp = './tests/test_cases/eofs_MJO.nc'
    check_corr_ds = xr.open_dataset(fp)
    check_corr_ds = whtools.interpolate_obs(check_corr_ds, MJO_obs.forecast_lons)

    check_modes = ['eof1_olr', 'eof2_olr', 'eof1_u200', 'eof2_u200', 'eof1_u850', 'eof2_u850']

    for cm in check_modes:
        ref = check_corr_ds[cm].values
        got = np.array(MJO_obs.MJO_fobs[cm])
        corr = np.corrcoef(ref, got)[0, 1]
        std_ratio = np.std(got) / np.std(ref)
        assert corr > 0.95, f"{cm}: correlation {corr:.3f} <= 0.95 (sign/structure error)"
        assert 0.85 < std_ratio < 1.18, f"{cm}: std ratio {std_ratio:.3f} (scaling error)"
        assert np.allclose(got, ref, atol=0.03), f"{cm}: values differ from reference beyond 0.03"

    return True

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

