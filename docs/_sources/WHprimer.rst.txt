Wheeler and Hendon Primer
==========================


A primer on the Wheeler and Hendon Real-time Multivariate MJO (RMM) Index
--------------------------------------------------------------------------

The Wheeler and Hendon Real-time Multivariate MJO (RMM) index is a widely used diagnostic tool for tracking and quantifying the Madden-Julian Oscillation (MJO), a significant tropical atmospheric phenomenon characterized by its eastward-propagating convective and atmospheric circulation anomalies.

Wheeler and Hendon Real-time Multivariate MJO (RMM) Index
----------------------------------------------------------

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