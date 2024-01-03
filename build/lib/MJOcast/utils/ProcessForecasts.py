# Import required libraries
import yaml  # For working with YAML configuration files
import os  # For operating system functionalities
import glob  # For file path pattern matching
import xarray as xr  # For working with labeled multi-dimensional arrays
import numpy as np  # For numerical operations
import eofs.standard as Eof_st  # For Empirical Orthogonal Function analysis
from eofs.multivariate.standard import MultivariateEof  # For multivariate EOF analysis
import matplotlib.pyplot as plt  # For creating plots
import matplotlib as mpl  # Matplotlib configuration and settings
import pandas as pd  # For working with data in tabular form
import copy  # For creating deep copies of objects
import re  # For regular expressions
from datetime import datetime  # For working with date and time
import sys  # For interacting with the Python interpreter
from MJOcast.utils.WHtools import (
    interpolate_obs,  # Custom utility functions for data manipulation
    check_or_create_paths,  # Custom utility functions for file and directory management
    convert_dates_to_string,  # Custom utility functions for date conversion
    check_lat_lon_coords,  # Custom utility functions for coordinate validation
    flip_lat_if_necessary,  # Custom utility functions for coordinate adjustment
    switch_lon_to_0_360,  # Custom utility functions for coordinate adjustment
    plot_phase_space  # Custom utility function for plotting MJO phase space
)

def make_DF_ense(files):
    """
    Create a DataFrame with files and their corresponding initialization dates.

    Parameters:
        files (list): List of file paths.

    Returns:
        DF (pd.DataFrame): DataFrame containing files and their initialization dates.
    """
    yr = []
    DF = pd.DataFrame({'File': files})

    # Initialize the 'Init' column with file names
    DF['Init'] = DF['File']

    for ee, File in enumerate(DF['File']):
        # Extract the initialization date from the file name using the provided function
        inst, matches = convert_dates_to_string(File)
        DF['Init'][ee] = matches[0]

    return DF


class MJOforecaster:
    """
    A class for forecasting the Madden-Julian Oscillation (MJO) using given parameters.

    This class provides methods for forecasting the MJO based on user-defined configurations.
    
    Parameters:
        yaml_file_path (str): The path to the YAML configuration file.
        eof_dict (dict): A dictionary containing Empirical Orthogonal Function (EOF) configurations.
        MJO_fobs (xarray.Dataset): Observed MJO data as an xarray Dataset.

    Attributes:
        yml_data (dict): Parsed data from the YAML configuration file.
        yml_usr_info (dict): User-defined information from the YAML configuration.
        forecast_lons (list): List of longitudes used for forecasting.
        base_dir (str): The base directory specified in the YAML configuration.
        eof_dict (dict): Dictionary containing EOF configurations.
        MJO_fobs (xarray.Dataset): Observed MJO data.
        made_forecast_file (bool): Indicates if a forecast file has been created.
    """
    def __init__(self,yaml_file_path,eof_dict,MJO_fobs):
        
        """
        Initialize the MJOforecaster instance with provided configurations.
        
        Parameters:
            yaml_file_path (str): The path to the YAML configuration file.
            eof_dict (dict): A dictionary containing Empirical Orthogonal Function (EOF) configurations.
            MJO_fobs (xarray.Dataset): Observed MJO data as an xarray Dataset.
        """
        
        with open(yaml_file_path, 'r') as file:
            yml_data = yaml.safe_load(file)
        yml_data
        
        self.yml_data = yml_data
        self.yml_usr_info = yml_data['user_defined_info']
        DSforexample = check_or_create_paths(yml_data)
        self.forecast_lons = DSforexample['lon']
        self.base_dir = yml_data['base_dir']
        self.eof_dict = eof_dict
        self.MJO_fobs = MJO_fobs
        
        
        self.made_forecast_file = False
        pass
    
    
    
    def get_forecast_LT_climo(self, yml_data, lons_forecast):
        """
        Get the forecast climatology dataset.

        Parameters:
            yml_data (dict): A dictionary containing user-defined information.
            lons_forecast (array-like): Array of forecast longitudes.

        Returns:
            DS_climo_forecast (xr.Dataset): Forecast climatology dataset.
        """
        # Extract the user-defined information from the input dictionary
        yml_usr_info = yml_data['user_defined_info']

        # Check if the user wants to use the forecast-dependent climatology
        if yml_usr_info['use_forecast_climo']:
            print('Using the forecast dependent climatology. Make sure you have generated it using ./Preprocessing_Scripts/*.ipynb.')
            # Load the forecast climatology dataset
            DS_climo_forecast = xr.open_dataset(self.base_dir+'/Forecast_Climo/Forecast_Climo.nc')
            #TODO add an option to generate/use their own forecast climo...
            DS_climo_forecast = interpolate_obs(DS_climo_forecast, lons_forecast)  # Check to make sure this works....
        else: 
            print('Using the climatology calculated by ERA5. It will be less skillful than a lead time dependent climatology.')
            print('Generate a lead time dependent climatology in ./Preprocessing_Scripts/*.ipynb for better results.')
            # Load the ERA5-based forecast climatology dataset
            obs_fp_Ec = os.path.join(self.base_dir,yml_usr_info['obs_data_loc'] + '/ERA5_climo.nc')
            DS_climo_forecast = xr.open_dataset(obs_fp_Ec)
            # Interpolate the ERA5-based climatology to match the forecast longitudes
            DS_climo_forecast = interpolate_obs(DS_climo_forecast, lons_forecast)  # Check to make sure this works....
        
        self.DS_climo_forecast = DS_climo_forecast
        return DS_climo_forecast


    def anomaly_LTD(self,yml_data, DS_CESM_for, DS_climo_forecast, numdays_out):
        """
        Calculate anomalies of u850, u200, and OLR from forecast and forecast climatology datasets.

        Parameters:
            yml_data (dict): A dictionary containing user-defined information.
            DS_CESM_for (xr.Dataset): Forecast dataset (CESM format).
            DS_climo_forecast (xr.Dataset): Forecast climatology dataset.
            numdays_out (int): Number of forecast days.

        Returns:
            U850_cesm_anom (xr.DataArray): Anomalies of u850.
            U200_cesm_anom (xr.DataArray): Anomalies of u200.
            OLR_cesm_anom (xr.DataArray): Anomalies of OLR.
        """
        yml_usr_info = yml_data['user_defined_info']

        u200vSTR = yml_usr_info['forecast_u200_name']
        u850vSTR = yml_usr_info['forecast_u850_name']
        olrvSTR = yml_usr_info['forecast_olr_name']

        #get the day of year of the forecasts
        fordoy = np.array(DS_CESM_for['time.dayofyear'])

        # Convert the 'fordoy' numpy array to a DataArray
        fordoy_da = xr.DataArray(fordoy, dims='lead', coords={'lead': range(len(fordoy))})
        # Use .sel() and .drop() with the DataArray fordoy
        DSclimo_doy = DS_climo_forecast.sel(dayofyear=fordoy_da,lead=slice(0,numdays_out-1))

        # u850:
        U850_cesm = DS_CESM_for[u850vSTR]
        U850_cesm_anom = xr.zeros_like(U850_cesm)
        temp_clim_u850 = np.array(DSclimo_doy['ua_850'])
        temp_clim_u850 = np.expand_dims(temp_clim_u850, 0)
        if temp_clim_u850.shape[1] == numdays_out + 1:
            temp_clim_u850 = temp_clim_u850[:, :numdays_out, :]
        U850_cesm_anom[:, :, :] = np.array(U850_cesm) - temp_clim_u850
        # u200:
        U200_cesm = DS_CESM_for[u200vSTR]
        U200_cesm_anom = xr.zeros_like(U200_cesm)
        temp_clim_u200 = np.array(DSclimo_doy['ua_200'])
        temp_clim_u200 = np.expand_dims(temp_clim_u200, 0)
        if temp_clim_u200.shape[1] == numdays_out + 1:
            temp_clim_u200 = temp_clim_u200[:, :numdays_out, :]
        U200_cesm_anom[:, :, :] = np.array(U200_cesm) - temp_clim_u200

        # OLR:
        OLRxr = DS_CESM_for[olrvSTR]
        OLR_cesm_anom = xr.zeros_like(DS_CESM_for[olrvSTR])
        temp_clim_olr = np.array(DSclimo_doy['rlut'])
        temp_clim_olr = np.expand_dims(temp_clim_olr, 0)
        if temp_clim_olr.shape[1] == numdays_out + 1:
            temp_clim_olr = temp_clim_olr[:, :numdays_out, :]
        OLR_cesm_anom[:, :, :] = np.array(OLRxr) - temp_clim_olr

        return U850_cesm_anom, U200_cesm_anom, OLR_cesm_anom


    def anomaly_ERA5(self,yml_data, DS_CESM_for, DS_climo_forecast, numdays_out):
        """
        Calculate anomalies of u850, u200, and OLR from forecast and ERA5 climatology datasets.

        Parameters:
            yml_data (dict): A dictionary containing user-defined information.
            DS_CESM_for (xr.Dataset): Forecast dataset (CESM format).
            DS_climo_forecast (xr.Dataset): ERA5 climatology dataset.
            numdays_out (int): Number of forecast days.

        Returns:
            U850_cesm_anom (xr.DataArray): Anomalies of u850.
            U200_cesm_anom (xr.DataArray): Anomalies of u200.
            OLR_cesm_anom (xr.DataArray): Anomalies of OLR.
        """
        yml_usr_info = yml_data['user_defined_info']

        u200vSTR = yml_usr_info['forecast_u200_name']
        u850vSTR = yml_usr_info['forecast_u850_name']
        olrvSTR = yml_usr_info['forecast_olr_name']
        
        obs_fp_Ec = os.path.join(self.base_dir,yml_usr_info['obs_data_loc'] + '/ERA5_climo.nc')
        ERA5clim = xr.open_dataset(obs_fp_Ec)
        U850_clim=ERA5clim['uwnd850'].to_dataset()
        U200_clim=ERA5clim['uwnd200'].to_dataset()
        OLR_clim=ERA5clim['olr'].to_dataset()

        #get the day of year of the forecasts
        fordoy = np.array(DS_CESM_for['time.dayofyear'])

        if fordoy[-1]>fordoy[0]:
            ### OLR ####
            OLRxr = DS_CESM_for[olrvSTR]
            OLR_cesm_anom = xr.zeros_like(DS_CESM_for[olrvSTR])
            temp_clim_olr = np.expand_dims(np.array(OLR_clim.sel(dayofyear=slice(fordoy[0],fordoy[-1]))['olr']),0)
            OLR_cesm_anom[:,:,:] = np.array(OLRxr)-temp_clim_olr

            ### u200 winds ####
            U200_cesm = DS_CESM_for[u200vSTR]
            U200_cesm_anom = xr.zeros_like(U200_cesm)
            temp_clim_u200 = np.expand_dims(np.array(U200_clim.sel(dayofyear=slice(fordoy[0],fordoy[-1]))['uwnd200']),0)
            U200_cesm_anom[:,:,:] = np.array(U200_cesm)-temp_clim_u200

            ### u850 winds ####
            U850_cesm = DS_CESM_for[u850vSTR]
            U850_cesm_anom = xr.zeros_like(U850_cesm)
            temp_clim_u850 = np.expand_dims(np.array(U850_clim.sel(dayofyear=slice(fordoy[0],fordoy[-1]))['uwnd850']),0)
            U850_cesm_anom[:,:,:] = np.array(U850_cesm)-temp_clim_u850

        else:
            print('...we crossed Jan 1...')
            ### OLR ####
            OLRxr = DS_CESM_for[olrvSTR]
            OLR_cesm_anom = xr.zeros_like(DS_CESM_for[olrvSTR])
            temp_clim_olr = np.concatenate([np.array(OLR_clim.sel(dayofyear=slice(fordoy[0],365))['olr']),np.array(OLR_clim.sel(dayofyear=slice(1,fordoy[-1]+1))['olr'])],axis=0)
            temp_clim_olr = np.expand_dims(temp_clim_olr,0)
            if temp_clim_olr.shape[1]==numdays_out + 1:
                temp_clim_olr = temp_clim_olr[:,:numdays_out,:,:]
            OLR_cesm_anom[:,:,:] = np.array(OLRxr)-temp_clim_olr

            ### u200 winds ####
            U200_cesm = DS_CESM_for[u200vSTR]
            U200_cesm_anom = xr.zeros_like(U200_cesm)
            temp_clim_u200 = np.concatenate([np.array(U200_clim.sel(dayofyear=slice(fordoy[0],365))['uwnd200']),np.array(U200_clim.sel(dayofyear=slice(1,fordoy[-1]+1))['uwnd200'])],axis=0)
            temp_clim_u200 = np.expand_dims(temp_clim_u200,0)
            if temp_clim_u200.shape[1]==numdays_out + 1:
                temp_clim_u200 = temp_clim_u200[:,:numdays_out,:,:]
            U200_cesm_anom[:,:,:] = np.array(U200_cesm)-temp_clim_u200

            ### u850 winds ####
            U850_cesm = DS_CESM_for[u850vSTR]
            U850_cesm_anom = xr.zeros_like(U850_cesm)
            temp_clim_u850 = np.concatenate([np.array(U850_clim.sel(dayofyear=slice(fordoy[0],365))['uwnd850']),np.array(U850_clim.sel(dayofyear=slice(1,fordoy[-1]+1))['uwnd850'])],axis=0)
            temp_clim_u850 = np.expand_dims(temp_clim_u850,0)
            if temp_clim_u850.shape[1]==numdays_out + 1:
                temp_clim_u850 = temp_clim_u850[:,:numdays_out,:,:]
            U850_cesm_anom[:,:,:] = np.array(U850_cesm)-temp_clim_u850

        return U850_cesm_anom, U200_cesm_anom, OLR_cesm_anom



    def filt_ndays(self,yml_data,DS_CESM_for,U850_cesm_anom,U200_cesm_anom,OLR_cesm_anom,DS_climo_forecast,numdays_out,AvgdayN,nensembs):
        """
        Perform anomaly filtering for atmospheric variables.

        Parameters:
            yml_data (dict): YAML data containing user-defined information.
            DS_CESM_for: Not defined in the code snippet provided.
            U850_cesm_anom (xarray.Dataset): Anomaly dataset for U850 variable.
            U200_cesm_anom (xarray.Dataset): Anomaly dataset for U200 variable.
            OLR_cesm_anom (xarray.Dataset): Anomaly dataset for OLR variable.
            DS_climo_forecast: Not defined in the code snippet provided.
            numdays_out: Not defined in the code snippet provided.
            AvgdayN: The number of days to be averaged for filtering.
            nensembs (int): Number of ensemble members.

        Returns:
            Updated U850_cesm_anom_filterd, U200_cesm_anom_filterd, OLR_cesm_anom_filterd datasets.
        """
        yml_usr_info = yml_data['user_defined_info']

        # Get variable names from user settings
        u200v = yml_usr_info['forecast_u200_name']
        u850v = yml_usr_info['forecast_u850_name']
        olrv = yml_usr_info['forecast_olr_name']

        # Define 120 days prior.
        first_date = U850_cesm_anom.time.values[0]- np.timedelta64(1, 'D') 
        first_date_120 =  U850_cesm_anom.time.values[0] - np.timedelta64(AvgdayN, 'D') 

        # Initialize arrays for filtered anomalies
        U850_cesm_anom_filterd = xr.zeros_like(U850_cesm_anom)
        U200_cesm_anom_filterd = xr.zeros_like(U200_cesm_anom)
        OLR_cesm_anom_filterd = xr.zeros_like(OLR_cesm_anom)

        ##get obs anomaly...
        obs_fp_MA = os.path.join(self.base_dir,yml_usr_info['obs_data_loc'] + '/ERA5_Meridional_Mean_Anomaly.nc')
        Obsanom = xr.open_dataset(obs_fp_MA)
        Obsanom = interpolate_obs(Obsanom, DS_CESM_for['lon'])
        OLR_anom = Obsanom['olr'].to_dataset().rename({'olr':yml_usr_info['forecast_olr_name']})
        U200_anom = Obsanom['uwnd200'].to_dataset().rename({'uwnd200':yml_usr_info['forecast_u200_name']})
        U850_anom = Obsanom['uwnd850'].to_dataset().rename({'uwnd850':yml_usr_info['forecast_u850_name']})

        for enen in range(nensembs):
            ### OLR anomaly filtering:
            tmpREolr=OLR_anom.sel(time=slice(first_date_120,first_date))
            tmpREolr=tmpREolr.drop_vars('dayofyear')
            fused_RE_for_OLR = xr.concat([tmpREolr,OLR_cesm_anom.sel(ensemble=enen).to_dataset()],dim='time')
            fused_RE_for_OLR_rolled = fused_RE_for_OLR.rolling(time=120, center=False,min_periods=1).mean().sel(time=slice(OLR_cesm_anom.time.values[0],OLR_cesm_anom.time.values[-1]))
            OLR_cesm_anom_filterd[enen,:,:]=OLR_cesm_anom.sel(ensemble=enen).values - fused_RE_for_OLR_rolled[olrv].values

            ### U200 anomaly filtering:
            tmpRE200=U200_anom.sel(time=slice(first_date_120,first_date))
            tmpRE200=tmpRE200.drop_vars('dayofyear')
            fused_RE_for_200 = xr.concat([tmpRE200,U200_cesm_anom.sel(ensemble=enen).to_dataset()],dim='time')
            fused_RE_for_200_rolled = fused_RE_for_200.rolling(time=120, center=False,min_periods=1).mean().sel(time=slice(U200_cesm_anom.time.values[0],U200_cesm_anom.time.values[-1]))
            U200_cesm_anom_filterd[enen,:,:]=U200_cesm_anom.sel(ensemble=enen).values - fused_RE_for_200_rolled[u200v].values

            ### U850 anomaly filtering:
            tmpRE850=U850_anom.sel(time=slice(first_date_120,first_date))
            tmpRE850=tmpRE850.drop_vars('dayofyear')
            fused_RE_for_850 = xr.concat([tmpRE850,U850_cesm_anom.sel(ensemble=enen).to_dataset()],dim='time')
            fused_RE_for_850_rolled = fused_RE_for_850.rolling(time=120, center=False,min_periods=1).mean().sel(time=slice(U850_cesm_anom.time.values[0],U850_cesm_anom.time.values[-1]))
            U850_cesm_anom_filterd[enen,:,:]=U850_cesm_anom.sel(ensemble=enen).values - fused_RE_for_850_rolled[u850v].values

        return OLR_cesm_anom_filterd,U200_cesm_anom_filterd,U850_cesm_anom_filterd


    def project_eofs(self,OLR_cesm_anom_filterd, U850_cesm_anom_filterd, U200_cesm_anom_filterd, numdays_out, nensembs, neofs_save, neof, eof_dict,svname,U200_cesm_anom):
        """
        Calculate and save RMM indices and EOFs.

        Parameters:
            OLR_cesm_anom_filterd (xarray.Dataset): Anomaly dataset for OLR variable.
            U850_cesm_anom_filterd (xarray.Dataset): Anomaly dataset for U850 variable.
            U200_cesm_anom_filterd (xarray.Dataset): Anomaly dataset for U200 variable.
            numdays_out (int): Number of days to project.
            nensembs (int): Number of ensemble members.
            neofs_save (int): Number of EOFs to save.
            neof (int): Number of EOFs to use for RMM calculation.
            eof_dict (dict): Dictionary containing normalization factors and other parameters.

        Returns:
            RMM1 (numpy.ndarray): RMM index 1.
            RMM2 (numpy.ndarray): RMM index 2.
            eofs_save (numpy.ndarray): Array of EOFs.
            sv_olr (numpy.ndarray): Scaled and normalized OLR data.
            sv_u200 (numpy.ndarray): Scaled and normalized U200 data.
            sv_u850 (numpy.ndarray): Scaled and normalized U850 data.
            sv_olr_unscaled (numpy.ndarray): Unscaled OLR data.
        """

        # Unpack the dictionary containing normalization factors and other parameters
        solver = eof_dict['solver']
        u200_norm = eof_dict['u200_norm']
        u850_norm = eof_dict['u850_norm']
        olr_norm = eof_dict['olr_norm']
        loc1 = eof_dict['loc1']
        loc2 = eof_dict['loc2']
        scale1 = eof_dict['scale1']
        scale2 = eof_dict['scale2']

        U850_cesm_anom_filterd_latmean = U850_cesm_anom_filterd
        U200_cesm_anom_filterd_latmean = U200_cesm_anom_filterd
        OLR_cesm_anom_filterd_latmean = OLR_cesm_anom_filterd

        # Initialize arrays for saving out data
        RMM1 = np.zeros([numdays_out, nensembs])
        RMM2 = np.zeros([numdays_out, nensembs])
        eofs_save = np.zeros([numdays_out, nensembs, neofs_save])
        sv_olr = np.zeros([numdays_out, nensembs, len(OLR_cesm_anom_filterd['lon'])])
        sv_u200 = np.zeros([numdays_out, nensembs, len(OLR_cesm_anom_filterd['lon'])])
        sv_u850 = np.zeros([numdays_out, nensembs, len(OLR_cesm_anom_filterd['lon'])])
        sv_olr_unscaled = np.zeros([numdays_out, nensembs, len(OLR_cesm_anom_filterd['lon'])])

        for enen in range(nensembs):
            # Normalize the anomaly data
            forc_u200_norm = np.array(U200_cesm_anom_filterd_latmean.sel(ensemble=enen) / u200_norm)
            forc_u850_norm = np.array(U850_cesm_anom_filterd_latmean.sel(ensemble=enen) / u850_norm)
            forc_olr_norm = np.array(OLR_cesm_anom_filterd_latmean.sel(ensemble=enen) / olr_norm)
            forc_olr_norm[:, -1] = 0
            forc_u200_norm[:, -1] = 0
            forc_u850_norm[:, -1] = 0

            neof = 2
            neofs_save = 15
            pj_sub = solver.projectField([forc_olr_norm, forc_u850_norm, forc_u200_norm], neofs=neofs_save)
            pj_saver_normalized = pj_sub / np.sqrt(solver.eigenvalues()[0:neofs_save])
            pj_sub = pj_sub[:, 0:neof] / np.sqrt(solver.eigenvalues()[0:neof])
            RMM1[:, enen] = pj_sub[:, loc1] * (scale1)
            RMM2[:, enen] = pj_sub[:, loc2] * (scale2)  # I think this is right

            sv_olr[:, enen, :] = forc_olr_norm.squeeze()
            sv_u200[:, enen, :] = forc_u200_norm.squeeze()
            sv_u850[:, enen, :] = forc_u850_norm.squeeze()
            eofs_save[:, enen, :] = pj_saver_normalized

        ###ensemble mean####
        forc_u200_norm = np.array((U200_cesm_anom_filterd_latmean/u200_norm).mean('ensemble'))
        forc_u850_norm = np.array((U850_cesm_anom_filterd_latmean/u850_norm).mean('ensemble'))
        forc_olr_norm = np.array((OLR_cesm_anom_filterd_latmean/olr_norm).mean('ensemble'))   
        forc_olr_norm[:,-1]=0
        forc_u200_norm[:,-1]=0
        forc_u850_norm[:,-1]=0
        pj_sub=solver.projectField([forc_olr_norm, forc_u850_norm, forc_u200_norm],neofs=neofs_save)
        pj_saver_normalized = pj_sub/np.sqrt(solver.eigenvalues()[0:neofs_save])
        pj_sub = pj_sub[:,0:neof]/np.sqrt(solver.eigenvalues()[0:neof])
        RMM1_emean = pj_sub[:,loc1]*(scale1) 
        RMM2_emean = pj_sub[:,loc2]*(scale2)  ### I think this is right
        ###ensemble mean####

        #grab the obs fields:
        obs_olr_norm=np.array(self.MJO_fobs.olr_norm.sel(time=slice(OLR_cesm_anom_filterd_latmean.time[0],OLR_cesm_anom_filterd_latmean.time[-1])))
        obs_u850_norm=np.array(self.MJO_fobs.u850_norm.sel(time=slice(OLR_cesm_anom_filterd_latmean.time[0],OLR_cesm_anom_filterd_latmean.time[-1])))
        obs_u200_norm=np.array(self.MJO_fobs.u200_norm.sel(time=slice(OLR_cesm_anom_filterd_latmean.time[0],OLR_cesm_anom_filterd_latmean.time[-1])))

        pj_sub_obs=solver.projectField([obs_olr_norm, obs_u850_norm, obs_u200_norm],neofs=neof)
        pj_sub_obs=pj_sub_obs/np.sqrt(solver.eigenvalues()[0:2])

        RMM1_obs_cera20c = pj_sub_obs[:,loc1]*(scale1) 
        RMM2_obs_cera20c = pj_sub_obs[:,loc2]*(scale2) 

        self.save_out_forecast_nc(RMM1,RMM2,RMM1_emean,RMM2_emean,RMM1_obs_cera20c,RMM2_obs_cera20c,eofs_save,self.MJO_fobs
                             ,sv_olr,sv_u200,sv_u850,eof_dict,neofs_save,OLR_cesm_anom_filterd_latmean,svname,U200_cesm_anom,U200_cesm_anom_filterd)

        return True

    def save_out_forecast_nc(self,RMM1,RMM2,RMM1_emean,RMM2_emean,RMM1_obs_cera20c,RMM2_obs_cera20c,eofs_save,MJO_fobs,
                             sv_olr,sv_u200,sv_u850,eof_dict,neofs_save,OLR_cesm_anom_filterd_latmean,svname,U200_cesm_anom,U200_cesm_anom_filterd):
        
        
        """
        Save forecasted MJO data to a netCDF file and set attribute information.

        This function saves forecasted MJO-related data into a netCDF file and assigns attribute
        information for better metadata representation.

        Parameters:
            RMM1 (numpy.ndarray): Array containing forecasted RMM1 data.
            RMM2 (numpy.ndarray): Array containing forecasted RMM2 data.
            RMM1_emean (numpy.ndarray): Array containing forecasted RMM1 ensemble mean data.
            RMM2_emean (numpy.ndarray): Array containing forecasted RMM2 ensemble mean data.
            RMM1_obs_cera20c (numpy.ndarray): Array containing observed RMM1 data (CERA-20C).
            RMM2_obs_cera20c (numpy.ndarray): Array containing observed RMM2 data (CERA-20C).
            eofs_save (numpy.ndarray): Array containing saved EOFs data.
            MJO_fobs (xarray.Dataset): Observed MJO data as an xarray Dataset.
            sv_olr (numpy.ndarray): Array containing saved normalized OLR data.
            sv_u200 (numpy.ndarray): Array containing saved normalized u200 data.
            sv_u850 (numpy.ndarray): Array containing saved normalized u850 data.
            eof_dict (dict): Dictionary containing EOF configurations.
            neofs_save (int): Number of EOFs to save.
            OLR_cesm_anom_filterd_latmean (xarray.DataArray): Filtered and latitude-mean OLR data.
            svname (str): Name of the netCDF file to save.
            U200_cesm_anom (numpy.ndarray): Array containing u200 anomalies data (CESM2).
            U200_cesm_anom_filterd (numpy.ndarray): Array containing filtered u200 anomalies data (CESM2).
        """

        solver=eof_dict['solver']

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


        MJO_for.to_netcdf(svname)
        print('saved: ',svname)  
        self.MJO_forecast_DS = MJO_for



    def check_forecast_files_runtime(self,for_file_list, yml_usr_info):
        """
        Check the forecast files for required variables and ensemble dimension.

        Parameters:
            for_file_list (str): File path of the forecast file to be checked.
            yml_usr_info (dict): Dictionary containing user settings from YAML file.

        Returns:
            Bingo (bool): True if the required variables are present, False otherwise.
            DS (xarray.Dataset): Updated dataset with added 'ensemble' dimension (if required).
        """
        convert_dates_to_string(for_file_list)
        # Open the forecast file as an xarray dataset
        DS = xr.open_dataset(for_file_list)

        #check for lat lon
        CheckMinLatLon = check_lat_lon_coords(DS)
        if not CheckMinLatLon:
            raise RuntimeError("the files MUST contain either a lat/lon or latitude/longitude coordinate") 

        #change variable name if necessary
        if 'longitude' in DS.coords:
            DS = DS.rename({'longitude':'lon','latitude':'lat'})

        #flip the orientation of the xarray if necessary 
        try:
            DS = flip_lat_if_necessary(DS)
            DS = switch_lon_to_0_360(DS)
        except:
            raise RuntimeError("it broke while re-orienting your forecast files to be S->N and 0-360 degrees lon") 

        # Get variable names from user settings
        u200v = yml_usr_info['forecast_u200_name']
        u850v = yml_usr_info['forecast_u850_name']
        olrv = yml_usr_info['forecast_olr_name']

        # Check if required variables are present in the dataset
        if u200v in DS.variables and u850v in DS.variables and olrv in DS.variables:
            Bingo = True
        else:
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            print('it looks like the defined u200,u850,or olr variables are not present in the file')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            Bingo = False

        # Get the ensemble dimension name from user settings
        ensemble_name = yml_usr_info['forecast_ensemble_dimension_name']

        # Check if the ensemble dimension is already present in the dataset
        if ensemble_name in DS.coords:
            print('ensemble dimension length:', len(DS[ensemble_name]))
            ense_length = len(DS[ensemble_name])
        else:
            # If the ensemble dimension is not present, try to add it
            try:
                # Create an 'ensemble' coordinate array and assign it to the dataset
                DS = DS.assign_coords({ensemble_name: np.arange(len(DS[ensemble_name]))})
                print('expanding coords to include ensemble')
                print('ensemble dimension length:', len(DS[ensemble_name]))
                ense_length = len(DS[ensemble_name])
            except:
                # If there is an error, add 'ensemble' as a new dimension to the dataset
                ensemble_values = [0]
                DS = DS.expand_dims(ensemble=ensemble_values)
                # Assign the coordinate values for the new 'ensemble' dimension
                DS = DS.assign_coords(ensemble=ensemble_values)
                ense_length = len(DS['ensemble'])

        if 'time' in DS.coords: 
            print(f"there are {len(DS['time.dayofyear'])} forecast lead days in these files")
            leaddays = len(DS['time.dayofyear'])
        else: 
            raise RuntimeError("the files MUST be compatible with xarray's DS['time.dayofyear'] functionality")

        return Bingo, DS, ense_length, leaddays


    def create_forecasts(self,num_files=None):
        """
        Create forecast files for a given range of latitude and settings.

        This function generates forecast files based on provided configurations, including latitude range,
        filtering parameters, and EOF settings. It processes each forecast file and saves the forecasted
        data into netCDF files. This function also performs data manipulation and filtering operations.

        Parameters:
            num_files (int, optional): The maximum number of forecast files to process. If None, process all files.

        Returns:
            DS_CESM_for (xr.Dataset): Forecast dataset containing processed forecasted data.
            OLR_cesm_anom_filterd (numpy.ndarray): Filtered OLR anomalies.
            U200_cesm_anom_filterd (numpy.ndarray): Filtered u200 anomalies.
            U850_cesm_anom_filterd (numpy.ndarray): Filtered u850 anomalies.
        """
        # Settings
        latwant = [16, -16]  # Latitudinal range
        AvgdayN = 120  # Filtering average
        neof = 2  # Number of EOFs (RMM indices)
        neofs_save = 15  # Number of total PCs to save

        # Extract information from the data dictionary
        yml_usr_info = self.yml_data['user_defined_info']
        datadir_Uwind = yml_usr_info['forecast_data_loc']

        # Get filenames in a dataframe
        FN_Uwind = sorted(glob.glob(yml_usr_info['forecast_data_loc'] + '/'+ yml_usr_info['forecast_data_name_str']))
        if len(FN_Uwind) == 0:
            raise FileNotFoundError(f"Files '{yml_usr_info['forecast_data_loc'] + '/'+ yml_usr_info['forecast_data_name_str']}' do not exist... "
                                    f"check your datestring of the filenames.")
            
            
        if num_files is None: 
            count_files = 100000000000000
        else: 
            count_files = 0

        #get the driver dataframe
        DF_Uwind = make_DF_ense(FN_Uwind)

        #get the climatology to create the anomalies
        DS_climo_forecast = self.get_forecast_LT_climo(self.yml_data, self.forecast_lons)

        #loop to make each forecast file:
        for FileCounter, eee in enumerate(range(0, len(DF_Uwind))):
            
            if count_files==num_files:
                print('done with requested number of files')
                break
            
            svname = yml_usr_info['output_files_loc'] + yml_usr_info['output_files_string'] + '_' + DF_Uwind['Init'][eee] + '.nc'

            if os.path.exists(svname):
                print(svname)
                print('The above forecast file already exists... Im skipping it... erase it if you want to make it again')
                continue  # Skip to the next iteration if the file exists

            #Load check and adjust the forecast file and get key variables
            print('reading forecasts file:', DF_Uwind['File'][eee])
            Bingo, DS_CESM_for, nensembs, numdays_out = self.check_forecast_files_runtime(DF_Uwind['File'][eee], yml_usr_info)   


            DS_CESM_for = DS_CESM_for.sel(lat =slice(latwant[1],latwant[0]))
            DS_CESM_for = DS_CESM_for.mean('lat')

            #initialize the files to save out:
            RMM1 = np.zeros([numdays_out, nensembs])
            RMM2 = np.zeros([numdays_out, nensembs])
            eofs_save = np.zeros([numdays_out, nensembs, neofs_save])
            sv_olr = np.zeros([numdays_out, nensembs, len(DS_CESM_for['lon'])])
            sv_u200 = np.zeros([numdays_out, nensembs, len(DS_CESM_for['lon'])])
            sv_u850 = np.zeros([numdays_out, nensembs, len(DS_CESM_for['lon'])])
            sv_olr_unscaled = np.zeros([numdays_out, nensembs, len(DS_CESM_for['lon'])])

            print('---- doing anomaly ----')
            try:
                if yml_usr_info['use_forecast_climo']:
                    print('im using the LT dp climo')
                    U850_cesm_anom,U200_cesm_anom,OLR_cesm_anom = self.anomaly_LTD(self.yml_data,DS_CESM_for,DS_climo_forecast,numdays_out)
                else: 
                    U850_cesm_anom,U200_cesm_anom,OLR_cesm_anom = self.anomaly_ERA5(self.yml_data,DS_CESM_for,DS_climo_forecast,numdays_out)
            except:
                raise RuntimeError("error happened while computing forecast runtime anomaly.. check the get_forecast_anom() function")
            print('---- done computing the anomaly----')


            print('--- filter out 120 days ---')
            try:
                OLR_cesm_anom_filterd,U200_cesm_anom_filterd,U850_cesm_anom_filterd=self.filt_ndays(self.yml_data,
                                                                                                    DS_CESM_for,
                                                                                                    U850_cesm_anom,
                                                                                                    U200_cesm_anom,
                                                                                                    OLR_cesm_anom,DS_climo_forecast,
                                                                                                    numdays_out,AvgdayN,nensembs)
                #function to run the anomaly... get_forecast_anom(yml_data,DS_CESM_for,DS_climo_forecast,MJO_for_obs)
            except:
                raise RuntimeError("error happened while computing filtering out the previous days.. check the filter_previous_days() function")
            print('--- done filtering out 120 days ---')


            print('--- project the EOFS ---')
            try:
                self.project_eofs(OLR_cesm_anom_filterd,
                                  U850_cesm_anom_filterd,
                                  U200_cesm_anom_filterd,
                                  numdays_out,nensembs,
                                  neofs_save,neof,
                                  self.eof_dict,
                                  svname,U200_cesm_anom)
                #function to run the anomaly... get_forecast_anom(yml_data,DS_CESM_for,DS_climo_forecast,MJO_for_obs)
            except:
                raise RuntimeError("error happened while computing filtering out the previous days.. check the filter_previous_days() function")
            print('--- done projecting the EOFS ---')

            self.DS_for = DS_CESM_for
            self.OLR_anom_filtered = OLR_cesm_anom_filterd
            self.U200_anom_filtered = U200_cesm_anom_filterd
            self.U850_anom_filtered = U850_cesm_anom_filterd
            self.made_forecast_file = True
            count_files+=1

        return DS_CESM_for,OLR_cesm_anom_filterd,U200_cesm_anom_filterd,U850_cesm_anom_filterd
    
    
    
    def plot_phase_space(self, Num_Ensembles, Lead):
        """
        Plot MJO RMM Ensemble phase space with a scatter plot indicating the progression of time.
        
        Parameters:
            date_start (str or pandas.Timestamp): Starting date for the plot.
            days_forward (int): Number of days to project forward.
        
        Raises:
            RuntimeError: If 'make_observed_MJO' hasn't been executed to create observational MJO data.
        """

        
        if not self.made_forecast_file:
            raise RuntimeError("create_forecasts MJO must have been created by running 'create_forecasts'")
            
            
        date_start = self.DS_for['time'][0]
        
        # Create a deep copy of the input date string and convert it to pandas.Timestamp
        date_start_str = str(date_start.values)[:10]
        date_start = pd.to_datetime(date_start_str)
        
        MJO_fobs = self.MJO_fobs
        DS_for = self.MJO_forecast_DS
        
        DS_for.sel(number=slice(0,Num_Ensembles))
        
        # Create a new figure and axis for the plot
        fig, ax = plt.subplots(figsize=(13, 10))
        fig.suptitle('MJO RMM phase space, Start Date: ' + date_start_str, fontsize=14)
        plt.title('')
        
        # Increment time by the specified number of days
        days_to_increment = Lead
        new_date = pd.to_datetime(date_start) + pd.Timedelta(days=Lead)
        
        # Define a colormap for the scatter plot
        cmap = plt.cm.plasma
        
        # Select RMM1 and RMM2 data for the specified date range
        pc1 = MJO_fobs['RMM1_obs'].sel(time=slice(date_start, new_date))
        pc2 = MJO_fobs['RMM2_obs'].sel(time=slice(date_start, new_date))
        
        time_values = np.arange(len(pc1))
        norm = mpl.colors.Normalize(vmin=0, vmax=len(time_values) - 1)
        
        # Plot the phase space diagram
        plot_phase_space(ax)
        
        # Create a scatter plot indicating time progression
        scatter = plt.scatter(pc1, pc2, c=time_values, cmap=cmap, s=50, norm=norm, alpha=0.75)
        colorbar = plt.colorbar(scatter, label='Time')
        colorbar.ax.tick_params(labelsize=14)  # Adjust fontsize of colorbar tick labels
        # Overlay black lines to connect scatter plot points
        plt.plot(pc1, pc2, color='black', alpha=0.4)
        
        
        pc1_for = DS_for['RMM1'].sel(time=slice(date_start, new_date))
        pc2_for = DS_for['RMM2'].sel(time=slice(date_start, new_date))
        plt.plot(pc1_for, pc2_for, color='red', alpha=0.3)
        pc1_em = DS_for['RMM1_emean'].sel(time=slice(date_start, new_date))
        pc2_em = DS_for['RMM2_emean'].sel(time=slice(date_start, new_date))
        plt.plot(pc1_em, pc2_em, color='red', alpha=0.8,linewidth=4)
        scatter = plt.scatter(pc1_em, pc2_em, color='red',s=50)
        
        # Save the plot to a file
        plt.savefig(self.yml_usr_info['output_plot_loc'] + '/' + './phase_space_MJO_forecast_' + date_start_str + '.png')
        
        # Close the plot
        plt.close()
