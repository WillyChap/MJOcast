# MJOcast
MJO Forecasting Repo

![example workflow](https://github.com/WillyChap/MJOcast/actions/workflows/pytest.yml/badge.svg)


## set up the docs script like so: 

https://www.youtube.com/watch?v=mV44dBi9qcQ


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
Thanks for the contributions from the community!


Requirements
------------
MJOcast is written for Python 3 (version >= 3.7) and depends on the packages NumPy, Pandas, SciPy, and Matplotlib. It runs on Linux
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
*Note for experienced users: We have slightly changed the API for the EOF calculation with version 1.4. to be more flexible 
for changes in the future. Please read the API docs or compare your code with the current example. The old API is still
working but will deprecate with one of the next releases. Adapting to the new interface will only take a few minutes.*

There are three basic entry points, of which you should read the documentation:

* Calculation of the EOFs: [calc_eofs_from_olr](https://willychap.github.io/MJOcast/api/omi_calculator.html#MJOcast.omi.omi_calculator.calc_eofs_from_olr).
* Calculation of the PCs: [calculate_pcs_from_olr](https://willychap.github.io/MJOcast/api/omi_calculator.html#MJOcast.omi.omi_calculator.calculate_pcs_from_olr).
* An OLR data container class, which has to be provided for the calculations: [OLRData](https://willychap.github.io/MJOcast/api/olr_handling.html#MJOcast.olr_handling.OLRData)

After you have installed MJOcast, you can download an
[example](https://github.com/willychap/MJOcast/tree/master/examples/) from the source, which consists of two files: 

* recalculate_original_omi.py: After downloading some data files, which are mentioned and linked in the source
  documentation of the example, you can run this example to recalculate the original OMI values. The script will save
  the computed Empirical Orthogonal Functions (EOFs) and the Principal Components (PCs) in two individual files, which
  can also be configured in the source code. In addition, it will save a few plots into a directory, which can
  also be configured in the source. These plots show the agreement with the original OMI values (slight deviations are 
  expected due to numerical differences. This is explained in detail in the corresponding software meta paper).

  Note that you can use this example also as a template to calculate OMI values with your own OLR data. 
  In order to do that, you only have to adapt two parts of the code, which are also marked in the code documentation.

  Note also that this script may run for one or two hours on common personal computer systems.

* evaluate_omi_reproduction.py: This script produces more detailed comparison plots and saves them into a directory.
  The script recalculate_original_omi.py has to be run first, since the evaluation script is based on the saved results.
  As for recalculate_original_omi.py, some file and directory names have to be adapted in the beginning of the code.

Both files are also available as Jupyter notebook files.

Documentation
-----------------
The documentation is found on [GitHub Pages](https://willychap.github.io/MJOcast/) and also in the docs
folder of the [source](https://github.com/willychap/MJOcast/docs/).

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
