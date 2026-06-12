# MJOcast
MJO Forecasting Repo

![example workflow](https://github.com/WillyChap/MJOcast/actions/workflows/pytest_obs.yml/badge.svg)


## set up the docs script like so: 

https://www.youtube.com/watch?v=mV44dBi9qcQ

Documentation
-----------------
The documentation is found on [GitHub Pages](https://willychap.github.io/MJOcast/) and also in the docs
folder of the [source](https://github.com/willychap/MJOcast/docs/).

Overview
--------

MJOcast: Empowering Atmospheric Scientists in MJO Subseasonal to Seasonal (S2S) Forecasting!

Welcome to MJOcast, a meticulously crafted Python package designed to support atmospheric scientists engaged in the captivating realm of Subseasonal to Seasonal (S2S) forecasting research.

MJOcast provides you with essential tools to effortlessly compute the Madden-Julian Oscillation (MJO) index, with a specific focus on the respected Wheeler and Hendon Real-Time Multivariate MJO Index (RMM). Additionally, it simplifies the generation of compiled netcdf variables and facilitates accurate forecasts of the MJO, complete with essential skill metrics. Particularly noteworthy is its ability to efficiently transform forecasted hindcast ensembles into MJO indices, making this often-complex task as straightforward as pointing the package at your model runs (with the appropriate model variables) and generating the ensemble of forecasts. Refer to the included examples to prepare your data. We've put considerable effort into ensuring compatibility with multiple common file formats from major modeling centers.

At present, we leverage the cutting-edge ERA5 reanalysis product for producing observations and performing the 120-day filtering of the forecast product. Each forecasted field is conveniently bundled with the forecast observations, enabling immediate skill validation. In hindcast mode, MJOcast capably handles forecasts generated from initializations spanning the years 1951 to 2023.

While the name MJOcast hints at future expansions to encompass a wider array of MJO indices, our initial release centers on the prominent RMM algorithm. We encourage you to explore the broader landscape of MJO toolboxes, including the noteworthy MJOindices toolbox. Although not tailored for forecasting, this alternative toolbox offers an array of other MJO index metrics.

For a comprehensive understanding of the scientific principles and methodologies that underpin MJOcast, we invite you to delve into the detailed exposition provided in Chapman et al. (2023). This publication forms the foundation for the concepts and capabilities that MJOcast brings to your research toolkit.

Citation
--------
If you use MJOcast in published research, please cite the corresponding paper: BLAH BLAH BLAH BLAH BLAH DOI: BLAH

Please check our [list of further scientific publications](https://willychap.github.io/MJOcast/references.html), on which the
implementation of the package is based. It is likely that some of these publications should also be cited.

Contributors
------------

We welcome any and all contributions to the communitty, please submit a github issue for increased functionalities! 

Requirements
------------
MJOcast is written for Python 3 (version >= 3.7) and depends on the packages xarray, eofs, NumPy, Pandas, SciPy, and Matplotlib (see setup and installation instructions). It runs on Linux
and Windows. Other operating systems have not been tested. 

There are optional requirements, which can be installed along with MJOcast (see below)

Installation
------------
MJOcast is available in the [Python Package Index (PyPI)](https://pypi.org/project/MJOcast/). It can be installed using, 
e.g., pip.
    
    pip3 install MJOcast
    
It can also be installed from the source, which is available on [Zenodo](http://dx.doi.org/10.5281/zenodo.3613752) and [GitHub](https://github.com/willychap/MJOcast). 
Download the source, move into the directory containing the file setup.py and run

    pip3 install .

The following optional dependency sets can additionaly be installed by adding ["set_name"] behind the above commands:
  * *full_func*: Install these packages to be sure that all options are really available. This might 
    require a higher Python version than for the core functionality alone. RMM may still be computed without these
    additional packages, but the number of alternative approaches is limited.
  * *dev*: packages that are needed for the further development, documentation and testing of MJOcast.

Examples: 

    pip3 install MJOcast[full_func]

or

    pip3 install -e .[dev]
 
Getting started / examples
--------------------------
MJOcast is driven by a YAML configuration file (see the `settings.yaml` examples in the source) and
has two stages. The forecast stage consumes the output of the observation stage:

```python
import MJOcast.utils.ProcessOBS as ProObs
import MJOcast.utils.ProcessForecasts as ProFo

# 1) Build the observed RMM EOFs/PCs (from ERA5 or your own obs file).
#    Writes MJO_obs_created.nc and verification plots into the configured output dirs.
MJO_obs = ProObs.MJOobsProcessor('settings.yaml')
OBS_DS, eof_list, pcs, MJO_fobs, eof_dict = MJO_obs.make_observed_MJO()

# 2) Project your forecast ensembles onto the observed EOFs to get RMM1/RMM2.
#    Writes one MJO forecast netCDF per initialization date.
MJO_for = ProFo.MJOforecaster('settings.yaml', MJO_obs.eof_dict, MJO_obs.MJO_fobs)
MJO_for.create_forecasts()
```

The main entry points are:

* `MJOcast.utils.ProcessOBS.MJOobsProcessor` and its `make_observed_MJO()` method — compute the
  observed multivariate EOFs and RMM indices. If the auto-derived EOF orientation does not correlate
  well with ERA5, pass a `scaling_dict` (`loc1`, `loc2`, `scale1`, `scale2`) to override it.
* `MJOcast.utils.ProcessForecasts.MJOforecaster` and its `create_forecasts()` method — compute
  lead-time anomalies, apply the 120-day filter, and project each ensemble member onto the observed
  EOFs.

**Note:** the 120-day filter requires the observed anomaly file
(`Observations/ERA5_Meridional_Mean_Anomaly.nc`) to cover the 120 days *before* each forecast
initialization date. If it does not (e.g. the bundled file ends before your init), `create_forecasts`
raises a clear error — regenerate/extend the file with `Preprocessing_Tools/Make_Obs.ipynb`.

The YAML file controls everything else: observation source (`use_era5`), the forecast variable names
as they appear in your netCDF files (`forecast_olr_name`, `forecast_u200_name`, `forecast_u850_name`),
the ensemble dimension name, the climatology mode, and all input/output paths. Runnable notebooks are
provided in `MJOcast/Example/` (with and without a pre-built YAML), and helper notebooks for preparing
inputs are in `MJOcast/Preprocessing_Tools/`.


Automated testing
-----------------
After you have installed MJOcast, you can also download
[unit and integration tests](https://github.com/cghoffmann/MJOcast/tree/master/tests/) from the source to check
your installation using pytest:

* Download the complete test directory to you local file system.

* Download the external input and reference data files from [Zenodo](https://doi.org/10.5281/zenodo.3746562). Details are given in a separate [Readme file](https://github.com/willychap/MJOcast/blob/master/tests/testdata/README). 
  * Note that some necessary data files are already included in the test directory in the repository. Make sure to download
    those files together with the tests. The data files on Zenodo are complementary and not 
    included in the repository for reasons of file size and ownership.

* Install the pytest package or simply install the MJOcast development dependencies with 

      pip3 install MJOcast[dev]

* Move into your local test directory and run

      pytest

* In the case that some tests are failing with FileNotFoundErrors, it is likely that the package code is actually working, but that the test 
  environment is not set up properly. You should check the following before contacting us:
   * Did you download the data files from the repository?
   * Did you download the data files from Zenodo?
   * Did you preserve the directory structure?
   * Did you execute pytest in the tests/ subdirectory (where pytest.ini is located)? 

* Note that the tests may run for a few hours on a common personal computer.
  * To get a quicker impression, you can omit slow tests by executing the following command. However, this will
    not check the core OMI computation, which is most important, of course.

        pytest -m 'not slow' 
