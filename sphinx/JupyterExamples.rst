Jupyter Notebook Examples
=========================

Link to Jupyter Notebook Examples: `Jupyter Examples <https://github.com/WillyChap/MJOcast/tree/main/MJOcast/Example>`_.

Start up
---------------
Follow these simple steps to quickly get started in a jupter notebook:

0. ** READ THIS FIRST:** :

This notebook shows the user how to drive the MJOcast on a single forecast file, and create desired output file. 

In this example, the output file has lead-time dependent bias removed, based on the lead time dependent climo file.

This forecast file climo file was made with teh Preprocessing_Tools, and MUST be completed prior to making this work on your files


1. **Install Packages:** :

.. code-block:: python

    import MJOcast.utils.ProcessForecasts as ProFo 
    import MJOcast.utils.ProcessOBS as ProObs
    import MJOcast.utils.WHtools as WHtools
    import importlib
    import xarray as xr
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib import image as mpimg
      
2. **Configuration:** Customize your forecasting experience by configuring the toolbox:

   2.1. *Manual Configuration:* Edit the YAML configuration file directly to define your forecast parameters, data paths, and other settings.

   2.2. *YAML Generator Tool:* Alternatively, use the provided YAML generator tool in the preprocessor Jupyter Notebook for a user-friendly configuration process.
   
   2.2 example:
   
.. code-block:: python
    
    # YAML file specifying output configuration
    output_yaml = 'output.yaml'

    # User information
    user_ = 'wchapman'

    # Base directory
    base_dir = '/glade/work/wchapman/MJOcast/MJOcast'

    # Flag to indicate whether to use ERA5 data
    use_era5 = True

    # Alternative name for user-defined observations (only applicable if user_era5=False)
    usr_named_obs = 'alternative_name_of_obs.nc'

    # Locations for observation and forecast data
    obs_data_loc = './Observations/'

    # Locations for observation and forecast data
    forecast_data_loc = '/glade/campaign/cgd/amp/wchapman/MJO_S2S_CESM2/p1/examples'

    # String pattern for MJO forecast data files
    forecast_data_name_str = 'S2Shindcast_cesm2cam6vs_MJOvars_*.nc'

    # Variable names for forecast data (and name of ensemble dimension)
    forecast_olr_name = 'rlut'
    forecast_u200_name = 'ua_200'
    forecast_u850_name = 'ua_850'
    forecast_ensemble_dimension_name = 'ensemble'

    # Location for diagnostic plots
    output_plot_loc = './output_plots/'

    # Location for output files
    output_files_loc = '/glade/campaign/cgd/amp/wchapman/MJO_S2S_CESM2/RMM_forecast/ERA5/'

    # Naming convention for MJO forecasts
    output_files_string = '/MJO_Forecast_Init'

    # Flag to indicate whether to use forecast climatology (turn False if you don't want to remove lead time dependent bias).
    use_forecast_climo = True

    # Flag to indicate whether to use observed climatology
    # (always false if "use_forecast_climo=true")
    use_observed_climo = False

    # Flag to regenerate climatology
    regenerate_climo = False

    # Flag to use Dask for climatology on a single node/machine
    use_dask_for_climo = True

    # Flag to regenerate climatology (single node/machine)
    regenerate_climo = False

    #create your driver yaml file: 
    WHtools.Create_Driver_Yaml('output.yaml', user_, base_dir, use_era5, usr_named_obs, obs_data_loc, forecast_data_loc, forecast_data_name_str, 
               forecast_olr_name, forecast_u200_name, forecast_u850_name, forecast_ensemble_dimension_name, 
               output_plot_loc, output_files_loc, output_files_string, use_forecast_climo, use_observed_climo, 
               regenerate_climo, use_dask_for_climo)


3. **Run your example and create your processed file: :** Begin forecasting:

.. code-block:: python

    #initialize the ObsProcessor:
    MJO_obs = ProObs.MJOobsProcessor(output_yaml)

    #Make the Observed MJO file:
    OBS_DS, eof_list, pcs, MJO_fobs, eof_dict = MJO_obs.make_observed_MJO()

    #Plot whatever day you would like for the obs phase space:
    MJO_obs.plot_phase_space('2001-01-01',60)

    #Now initialize the 
    MJO_for = ProFo.MJOforecaster(output_yaml,MJO_obs.eof_dict,MJO_obs.MJO_fobs)

    #Now create the forecasts
    DS_CESM_for,OLR_cesm_anom_filterd,U200_cesm_anom_filterd,U850_cesm_anom_filterd = MJO_for.create_forecasts(num_files=1)

    #plot the phase space diagram for your forecast: 
    MJO_for.plot_phase_space(12,15)
    Load the YAML configuration file:

.. code-block:: python

    [output]:
    Number of forecast files to process: 1
    expanding coords to include ensemble
    ensemble dimension length: 11
    there are 46 forecast lead days in these files
    Initial look at forecast files passes the first test
     ----- Taking the EOF ----- 
    ...done making observed EOFS, check ./output_plots/*.png for verification metrics...
    ...attaching the BOM index for verification...
    ...saved OLR obs file...
    Number of forecast files to process: 1
    expanding coords to include ensemble
    ensemble dimension length: 11
    there are 46 forecast lead days in these files
    Initial look at forecast files passes the first test
    Using the forecast dependent climatology. Make sure you have generated it using ./Preprocessing_Scripts/*.ipynb.
    reading forecasts file: /glade/campaign/cgd/amp/wchapman/MJO_S2S_CESM2/p1/examples/S2Shindcast_cesm2cam6vs_MJOvars_25jan2018.nc
    expanding coords to include ensemble
    ensemble dimension length: 11
    there are 46 forecast lead days in these files
    ---- doing anomaly ----
    im using the LT dp climo
    ---- done computing the anomaly----
    --- filter out 120 days ---
    --- done filtering out 120 days ---
    --- project the EOFS ---
    saved:  /glade/campaign/cgd/amp/wchapman/MJO_S2S_CESM2/RMM_forecast/ERA5//MJO_Forecast_Init_25jan2018.nc
    --- done projecting the EOFS ---

A created dataset is then saved here: 

.. code-block:: python

    DSdone = xr.open_dataset('/glade/campaign/cgd/amp/wchapman/MJO_S2S_CESM2/RMM_forecast/ERA5//MJO_Forecast_Init_25jan2018.nc')

Output figures are saved here: 

.. code-block:: python
    
    "./output_plots/observed_eofs.png"

Please visit the provided links at the top of this page to see the running jupyter notebooks! 

Follow these steps to make the most of MJOforecaster in your projects. Happy forecasting!