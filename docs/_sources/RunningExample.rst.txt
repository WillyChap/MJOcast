Run Examples
=============

Start up
---------------
Follow these simple steps to quickly get started:

1. **Installation:** Install the toolbox using the provided *installation instructions* or simply type the following command in your terminal:

.. code-block:: bash

    pip install MJOcast
      
2. **Configuration:** Customize your forecasting experience by configuring the toolbox:

   2.1. *Manual Configuration:* Edit the YAML configuration file directly to define your forecast parameters, data paths, and other settings.

   2.2. *YAML Generator Tool:* Alternatively, use the provided YAML generator tool in the preprocessor Jupyter Notebook for a user-friendly configuration process.

3. **Quick Start Usage:** Begin forecasting with ease:

   Import the necessary modules:

.. code-block:: python

    import MJOcast.utils.ProcessForecasts as ProFo 
    import MJOcast.utils.ProcessOBS as ProObs

Load the YAML configuration file:

.. code-block:: python

    yaml_file_path = './settings.yaml'

Create Observational MJO and EOFs [from either user specified data or provide ERA5]:

.. code-block:: python

    MJO_obs = ProObs.MJOobsProcessor(yaml_file_path)
    OBS_DS, eof_list, pcs, MJO_fobs, eof_dict = MJO_obs.make_observed_MJO()

Generate Hindcast .nc Files for each forecast:
    
.. code-block:: python

    MJO_for = ProFo.MJOforecaster(yaml_file_path, MJO_obs.eof_dict, MJO_obs.MJO_fobs)
    DS_CESM_for, OLR_cesm_anom_filtered, U200_cesm_anom_filtered, U850_cesm_anom_filtered = MJO_for.create_forecasts(num_files=1)

4. **Documentation:** For in-depth guidance, explore the code API documentation. It provides detailed usage instructions, function explanations, and examples.

Follow these steps to make the most of MJOforecaster in your projects. Happy forecasting!