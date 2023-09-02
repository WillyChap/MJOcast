.. MJOcast documentation master file, created by
   sphinx-quickstart on Tue Aug 15 19:31:55 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MJOcast's documentation!
===================================

MJOforecaster Documentation
===========================

Welcome to the official documentation for the MJOforecaster toolbox - your comprehensive python package for generating MJO hindcast forecasts based on Wheeler and Hendon 2004 Real-time Multivariate MJO (RMM) methodology. For a primer on the WH RMM index, and it's forecasting see - :doc:`WHprimer`.

Overview
--------

The MJOforecaster aims to simplify the creation of precise MJO hindcast forecasts. Leveraging empirical orthogonal function (EOF) analyses and specialized data processing techniques, this toolbox empowers users to generate MJO forecasts within user-defined parameters. This documentation offers in-depth guidance on the toolbox's functionalities, usage, and advanced features.

To gain a clear understanding of the process for generating hindcasts of the RMM index, refer to the "examples" Jupyter notebooks available in the **MJOcast** GitHub repository. These notebooks provide a practical showcase, guiding users through the steps to master the art of forecast creation.

Additionally, MJOcast provides the functionallity to make publication ready figures. See the *figures page* to see the provided figure options.

When the MJOcast toolbox is run an individual MJO RMM forecast file will be created for each provided forecast. This file will also contain the observed RMM index for quick forecast assessment. The generated .nc file will look like so: 

.. code-block:: python

    ncdump -h *.nc 
    
    MJO_for = xr.Dataset(
        {
            "RMM1": (["time","number"],RMM1),
            "RMM2": (["time","number"],RMM2),
            "RMM1_emean": (["time"],RMM1_emean),
            "RMM2_emean": (["time"],RMM2_emean),
            "RMM1_obs":(["time"],RMM1_obs_cera20c),
            "RMM2_obs":(["time"],RMM2_obs_cera20c),
            "eofs_save":(["time","number",'neigs'],eofs_save),
            "OLR_norm":(["time","number","longitude"],sv_olr), 
            "eof1_olr":(["longitude"],np.array(self.MJO_fobs['eof1_olr'])),
            "eof2_olr":(["longitude"],np.array(self.MJO_fobs['eof2_olr'])),
            "eof1_u850":(["longitude"],np.array(self.MJO_fobs['eof1_u850'])),
            "eof2_u850":(["longitude"],np.array(self.MJO_fobs['eof2_u850'])),
            "eof1_u200":(["longitude"],np.array(self.MJO_fobs['eof2_u200'])),
            "eof2_u200":(["longitude"],np.array(self.MJO_fobs['eof2_u200'])),
            "u200_norm":(["time","number","longitude"],sv_u200), 
            "u850_norm":(["time","number","longitude"],sv_u850),
            "U200_cesm_anom":(["number","time","longitude"],np.array(U200_cesm_anom)), 
            "U200_cesm_anom_filterd":(["number","time","longitude"],np.array(U200_cesm_anom_filterd)), 
            "eig_vals":(['neigs'],solver.eigenvalues(neigs=neofs_save))

        },
        coords={
            "time":OLR_cesm_anom_filterd_latmean.time,
            "longitude":np.array(OLR_cesm_anom_filterd_latmean.lon),
            "neigs":np.arange(0,neofs_save)
        },)


        MJO_for.attrs["title"] = "MJO RMM Forecast the projected eof(u850,u200,OLR)"
                MJO_for.attrs["description"] = "MJO Forecast in the Prescribed Forecast dataset calculated as in Wheeler and Hendon 2004, a 120-day filter, 15S-15N averaged variables."
        MJO_for.attrs["notes"] = "ONLY Variables RMM1 and RMM2 have been flipped and switched -from eofs_save- to match standard MJO conventions."

        MJO_for.RMM1.attrs['units'] = 'stddev'
        MJO_for.RMM1.attrs['standard_name'] = 'RMM1'
        MJO_for.RMM1.attrs['long_name'] = 'Real-time Multivariate MJO Index 1'

        MJO_for.RMM2.attrs['units'] = 'stddev'
        MJO_for.RMM2.attrs['standard_name'] = 'RMM2'
        MJO_for.RMM2.attrs['long_name'] = 'Real-time Multivariate MJO Index 2'

        MJO_for.RMM1_emean.attrs['units'] = 'stddev'
        MJO_for.RMM1_emean.attrs['long_name'] = 'ensemble mean RMM1 forecast'

        MJO_for.RMM2_emean.attrs['units'] = 'stddev'
        MJO_for.RMM2_emean.attrs['long_name'] = 'ensemble mean RMM2 forecast'

        MJO_for.RMM1_obs.attrs['units'] = 'stddev'
        MJO_for.RMM1_obs.attrs['long_name'] = 'Observed RMM1 (ERA-5 or provided obs file)'

        MJO_for.RMM2_obs.attrs['units'] = 'stddev'
        MJO_for.RMM2_obs.attrs['long_name'] = 'Observed RMM2 (ERA-5 or provided obs file)'

        MJO_for.eofs_save.attrs['long_name'] = 'Empirical Orthogonal Functions (EOFs) of MJO observed variables'
        MJO_for.eofs_save.attrs['description'] = 'Matrix containing the EOFs of MJO variables for each ensemble member'

        MJO_for.OLR_norm.attrs['units'] = 'stddev'
        MJO_for.OLR_norm.attrs['standard_name'] = 'forecast OLR - normalized'
        MJO_for.OLR_norm.attrs['long_name'] = 'Normalized Outgoing Longwave Radiation'

        MJO_for.u200_norm.attrs['units'] = 'stddev'
        MJO_for.u200_norm.attrs['standard_name'] = 'forecast u200 normalized'
        MJO_for.u200_norm.attrs['long_name'] = 'Normalized Zonal Wind at 200mb'

        MJO_for.u850_norm.attrs['units'] = 'stddev'
        MJO_for.u850_norm.attrs['standard_name'] = 'forecasts u850 normalized'
        MJO_for.u850_norm.attrs['long_name'] = 'Normalized Zonal Wind at 850mb'

        MJO_for.U200_cesm_anom.attrs['long_name'] = 'Zonal Wind Anomalies at 200mb'
        MJO_for.U200_cesm_anom_filterd.attrs['long_name'] = '120 day Filtered Anomalies of all variables'

        MJO_for.eig_vals.attrs['long_name'] = 'Eigenvalues of Empirical Orthogonal Functions (EOFs)'
        MJO_for.eig_vals.attrs['description'] = 'Eigenvalues corresponding to the EOFs of MJO variables'


Key Features
------------
- **Hindcast Forecasts:** Generate MJO hindcast forecasts using provided configurations.
- **Data Manipulation:** Process, filter, and manipulate input forecast data for accurate hindcasting.
- **Empirical Orthogonal Functions (EOFs):** Utilize EOF analysis to identify MJO patterns and generate forecasts.
- **Metadata Management:** Set attribute information for generated forecast files.
- **Comprehensive Customization:** Customize forecasts by adjusting parameters, latitudinal range, and more.

What You Need to Get Started
----------------------------

To begin your journey with MJOcast, make sure you have the following:

1. Hindcast Files: Organize your hindcast files into different initializations, and name them with the initialization time (e.g., `CESM2_init_May032002.nc`).
   
   1.1. Forecast Files: These should encompass all ensembles and can be arranged in various grid orientations (such as lat, lon, ensemble, time or time, ensemble, lat, lon, etc.).
   
   1.2. Required Variables: Forecast files must include zonal wind at 200mb, zonal wind at 850mb, and Outgoing Longwave Radiation.
   
   1.3. Geographical Orientation: Files can be either in the -180°E to 180°W range or 0° to 360°W, and oriented from 90°S to 90°N or 90°N to 90°S.

2. Positive Attitude: A can-do attitude will make your forecasting experience even more rewarding!

**Optional**

3. User-Provided Observations: If you have your own observations, you can integrate them.

   3.1. Observation Generation: Consult the provided **preprocessing folder** in the code repository to generate observations from your own .nc files.
   
   3.2. Default Option: Alternatively, utilize the provided ERA5 observations (spanning 1950-2022) from the GitHub repository.

4. User-Provided Lead-Time Dependent Climatology: If available, incorporate your own climatology.

   4.1. Climatology Generation: Explore the **preprocessing folder** in the code repository to create this lead-time dependent climatology.
   
   4.2. Value of Climatology: Discover the benefits of using this climatology as explained below. See - :doc:`WhyLead` for the details.
   
   4.3. Default Option: Alternatively, the ERA5 data (or provided observations) will be used.

Prepare these essentials to embark on your journey with MJOforecaster and unlock its full potential.


Getting Started
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

Contents
--------

The documentation is organized into several sections:

- **Installation:** Details on how to install the MJOforecaster toolbox.
- **Configuration:** Information on setting up the YAML configuration file.
- **Usage Guide:** A comprehensive guide on using the toolbox's functionalities.
- **Advanced Features:** Exploring advanced features and customization options.
- **API Reference:** Detailed information about classes, methods, and parameters.
- **Examples:** Real-world examples illustrating different aspects of the toolbox's usage.

A primer on the Wheeler and Hendon Real-time Multivariate MJO (RMM) Index
=========================================================================

The Wheeler and Hendon Real-time Multivariate MJO (RMM) index is a widely used diagnostic tool for tracking and quantifying the Madden-Julian Oscillation (MJO), a significant tropical atmospheric phenomenon characterized by its eastward-propagating convective and atmospheric circulation anomalies.

Wheeler and Hendon Real-time Multivariate MJO (RMM) Index
=========================================================

The Wheeler and Hendon Real-time Multivariate MJO (RMM) index is a widely used diagnostic tool for tracking and quantifying the Madden-Julian Oscillation (MJO), a significant tropical atmospheric phenomenon characterized by its eastward-propagating convective and atmospheric circulation anomalies.

Construction of the Wheeler and Hendon RMM Index
------------------------------------------------

The RMM index is created through a principal component analysis (PCA) of key meteorological variables that capture the MJO's signature across the tropical Indian Ocean and Western Pacific regions. These variables include:

1. **Outgoing Longwave Radiation (OLR):** A measure of Earth's thermal radiation emitted to space. Negative OLR anomalies indicate enhanced convection associated with the active phase of the MJO.

2. **Zonal Wind at 850mb and 200mb:** Zonal wind anomalies at these pressure levels indicate changes in atmospheric circulation due to the MJO. The phase relationship between these winds determines the MJO's movement.

The PCA results in two principal components (PCs), known as RMM1 and RMM2. RMM1 represents the MJO's amplitude, while RMM2 captures its phase. Prior to PCA, anomalies of these variables are calculated by removing the 120-day mean. This process helps eliminate low-frequency variability, such as the effects of El Niño-Southern Oscillation (ENSO).

The RMM index is formed by combining these PCs into a complex time series:

- **RMM Index = RMM1 + i * RMM2**

Interpretation of the RMM Index
-------------------------------

- **Amplitude (RMM1):** Larger RMM1 values indicate stronger MJO activity. Negative values correspond to the suppressed phase (inactive MJO), and positive values indicate the active phase (enhanced convection and circulation).

- **Phase (RMM2):** The RMM2 phase determines the geographic location of the MJO's convection center. It's categorized into eight phases, each linked to specific equatorial regions.

The Wheeler and Hendon RMM index is a valuable tool for monitoring and forecasting the MJO's behavior, aiding in predictions of tropical weather patterns, precipitation, and atmospheric conditions on various time scales.

Note: This information is based on knowledge available up to September 2021. For the latest updates or changes, consult authoritative sources.

Forecasting with the RMM System
-------------------------------

To effectively forecast with the RMM system, certain preprocessing steps are required to prepare both observational and forecast data. Here's an overview of the necessary procedures:

1. **Incorporating the Past 120 Days (Reanalysis):** 
   The reanalysis data, such as ERA5, plays a pivotal role in refining forecasts. It's crucial to append the preceding 120 days of reanalysis data to the forecast data. This historical context helps remove the 120-day mean, effectively eliminating low-frequency variability (e.g., ENSO effects) from the analysis.

2. **Constructing Anomaly Fields:**
   Once the forecast and appended reanalysis data are available, compute the anomaly fields by subtracting the climatological mean from each day's data. This process isolates the deviations from the typical atmospheric conditions, allowing a more accurate analysis of MJO-related variability.

3. **Projecting ERA5 Empirical Orthogonal Functions (EOFS) Modes:**
   The next step involves projecting the modes of the Empirical Orthogonal Functions (EOFS) obtained from the ERA5 dataset onto the anomaly fields of the forecast data. This projection aligns the MJO-related variability between the observational and forecast datasets, facilitating comparison and enhancing forecast accuracy.

By following these steps, the RMM system harmonizes historical and forecast data, removing spurious signals, and ensuring that the forecast captures the genuine MJO behavior. This preprocessing lays the foundation for robust and accurate MJO hindcast forecasts.

.. _Why_leadtime:

Why Provide a Lead-Time Dependent Model Climatology in S2S Forecasting?
-----------------------------------------------------------------------

In Subseasonal-to-Seasonal (S2S) forecasting, incorporating a lead-time dependent model climatology offers several advantages that contribute to improved forecast accuracy and reliability. Here's why it's considered a beneficial practice:

1. **Reduced Biases:** A standard climatology applied uniformly across all lead times might not account for the dynamic nature of atmospheric conditions, especially as a model slips into its own biased climatology, which is distinct from the real-world climatology. A lead-time dependent model climatology adapts to the evolving climate, mitigating biases that arise from using fixed observed climatological values.

2. **Capturing Evolution:** The climate system's behavior changes as forecasts extend further into the future. A lead-time dependent climatology better captures this evolving nature, ensuring that the forecast is aligned with the expected conditions at various lead times.

3. **Improved Skill:** Incorporating a climatology that considers the lead time enhances the forecast's skill by accounting for the removal of the biased evolution of the model's climate. This leads to better capturing the evolving patterns, resulting in more accurate predictions.

4. **Contextual Insights:** A lead-time dependent climatology provides contextual insights into how the model's climate system evolves over time. This understanding enhances the interpretation of forecast anomalies and aids in identifying significant departures from the norm.

In summary, integrating a lead-time dependent model climatology acknowledges the changing nature of the model's climate system and how it transitions into its own attractor space. By harnessing its variability, this approach enhances the accuracy and reliability of S2S forecasts. By adapting to evolving atmospheric conditions, this approach helps create more skillful and informative predictions.


Feedback and Support
---------------------

We value your feedback! If you have any questions, suggestions, or issues while using the MJOforecaster toolbox, don't hesitate to reach out to us via `GitHub Issues <https://github.com/WillyChap/MJOcast/issues>`_.

Get ready to revolutionize your MJO hindcasting process with the MJOforecaster toolbox!
This is a descriptions of this website... 

.. toctree::
   :maxdepth: 4
   :caption: Contents:
   
   gettingstarted
   WHprimer
   WhyLead
   RunningExample
   JupyterExamples
   summary
   code

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
