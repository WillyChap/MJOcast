Getting Started
================

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
   
   4.2. Value of Climatology: Discover the benefits of using this climatology as explained below. See :any:`Why_leadtime` for the details.
   
   4.3. Default Option: Alternatively, the ERA5 data (or provided observations) will be used.

Prepare these essentials to embark on your journey with MJOforecaster and unlock its full potential.
