import math
from pathlib import Path
import os
import pytest

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

def test_default_creation():
    '''Returns a MJO_obs instance'''
    assert create_MJOobs()

def test_default_creation():
    '''Returns a MJO_obs instance'''
    assert create_MJOobs_file()
    