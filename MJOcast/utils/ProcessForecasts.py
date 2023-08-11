import yaml
import os
import glob
import xarray as xr
import numpy as np
import eofs.standard as Eof_st
from eofs.multivariate.standard import MultivariateEof
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import copy
import re
from datetime import datetime
import os
import sys
from MJOcast.utils.WHtools import (
    interpolate_obs,
    check_or_create_paths,
    convert_dates_to_string,
    check_lat_lon_coords,
    flip_lat_if_necessary,
    switch_lon_to_0_360
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
    def __init__(self,yaml_file_path,eof_dict,MJO_fobs):
        
        #Global attributes that the functions need! 
        
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
            DS_climo_forecast = xr.open_dataset('./Forecast_Climo/Forecast_Climo.nc')
            #TODO add an option to generate/use their own forecast climo...
        else: 
            print('Using the climatology calculated by ERA5. It will be less skillful than a lead time dependent climatology.')
            print('Generate a lead time dependent climatology in ./Preprocessing_Scripts/*.ipynb for better results.')
            # Load the ERA5-based forecast climatology dataset
            DS_climo_forecast = xr.open_dataset('./Observations/ERA5_climo.nc')
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
        DSclimo_doy = DS_climo_forecast.sel(dayofyear=fordoy_da)

        # u850:
        U850_cesm = DS_CESM_for[u850vSTR].squeeze()
        U850_cesm_anom = xr.zeros_like(U850_cesm.squeeze())
        temp_clim_u850 = np.array(DSclimo_doy['ua_850'])
        temp_clim_u850 = np.expand_dims(temp_clim_u850, 0)
        if temp_clim_u850.shape[1] == numdays_out + 1:
            temp_clim_u850 = temp_clim_u850[:, :numdays_out, :]
        U850_cesm_anom[:, :, :] = np.array(U850_cesm) - temp_clim_u850

        # u200:
        U200_cesm = DS_CESM_for[u200vSTR].squeeze()
        U200_cesm_anom = xr.zeros_like(U200_cesm.squeeze())
        temp_clim_u200 = np.array(DSclimo_doy['ua_200'])
        temp_clim_u200 = np.expand_dims(temp_clim_u200, 0)
        if temp_clim_u200.shape[1] == numdays_out + 1:
            temp_clim_u200 = temp_clim_u200[:, :numdays_out, :]
        U200_cesm_anom[:, :, :] = np.array(U200_cesm) - temp_clim_u200

        # OLR:
        OLRxr = DS_CESM_for[olrvSTR].squeeze()
        OLR_cesm_anom = xr.zeros_like(DS_CESM_for[olrvSTR].squeeze())
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

        ERA5clim = xr.open_dataset('./Observations/ERA5_climo.nc')
        U850_clim=ERA5clim['uwnd850'].to_dataset()
        U200_clim=ERA5clim['uwnd200'].to_dataset()
        OLR_clim=ERA5clim['olr'].to_dataset()

        #get the day of year of the forecasts
        fordoy = np.array(DS_CESM_for['time.dayofyear'])

        if fordoy[-1]>fordoy[0]:
            ### OLR ####
            OLRxr = DS_CESM_for[olrvSTR].squeeze()
            OLR_cesm_anom = xr.zeros_like(DS_CESM_for[olrvSTR].squeeze())
            temp_clim_olr = np.expand_dims(np.array(OLR_clim.sel(dayofyear=slice(fordoy[0],fordoy[-1]))['olr']),0)
            OLR_cesm_anom[:,:,:] = np.array(OLRxr)-temp_clim_olr

            ### u200 winds ####
            U200_cesm = DS_CESM_for[u200vSTR].squeeze()
            U200_cesm_anom = xr.zeros_like(U200_cesm.squeeze())
            temp_clim_u200 = np.expand_dims(np.array(U200_clim.sel(dayofyear=slice(fordoy[0],fordoy[-1]))['uwnd200']),0)
            U200_cesm_anom[:,:,:] = np.array(U200_cesm)-temp_clim_u200

            ### u850 winds ####
            U850_cesm = DS_CESM_for[u850vSTR].squeeze()
            U850_cesm_anom = xr.zeros_like(U850_cesm.squeeze())
            temp_clim_u850 = np.expand_dims(np.array(U850_clim.sel(dayofyear=slice(fordoy[0],fordoy[-1]))['uwnd850']),0)
            U850_cesm_anom[:,:,:] = np.array(U850_cesm)-temp_clim_u850

        else:
            print('...we crossed Jan 1...')
            ### OLR ####
            OLRxr = DS_CESM_for[olrvSTR].squeeze()
            OLR_cesm_anom = xr.zeros_like(DS_CESM_for[olrvSTR].squeeze())
            temp_clim_olr = np.concatenate([np.array(OLR_clim.sel(dayofyear=slice(fordoy[0],365))['olr']),np.array(OLR_clim.sel(dayofyear=slice(1,fordoy[-1]+1))['olr'])],axis=0)
            temp_clim_olr = np.expand_dims(temp_clim_olr,0)
            if temp_clim_olr.shape[1]==numdays_out + 1:
                temp_clim_olr = temp_clim_olr[:,:numdays_out,:,:]
            OLR_cesm_anom[:,:,:] = np.array(OLRxr)-temp_clim_olr

            ### u200 winds ####
            U200_cesm = DS_CESM_for[u200vSTR].squeeze()
            U200_cesm_anom = xr.zeros_like(U200_cesm.squeeze())
            temp_clim_u200 = np.concatenate([np.array(U200_clim.sel(dayofyear=slice(fordoy[0],365))['uwnd200']),np.array(U200_clim.sel(dayofyear=slice(1,fordoy[-1]+1))['uwnd'])],axis=0)
            temp_clim_u200 = np.expand_dims(temp_clim_u200,0)
            if temp_clim_u200.shape[1]==numdays_out + 1:
                temp_clim_u200 = temp_clim_u200[:,:numdays_out,:,:]
            U200_cesm_anom[:,:,:] = np.array(U200_cesm)-temp_clim_u200

            ### u850 winds ####
            U850_cesm = DS_CESM_for[u850vSTR].squeeze()
            U850_cesm_anom = xr.zeros_like(U850_cesm.squeeze())
            temp_clim_u850 = np.concatenate([np.array(U850_clim.sel(dayofyear=slice(fordoy[0],365))['uwnd850']),np.array(U850_clim.sel(dayofyear=slice(1,fordoy[-1]+1))['uwnd'])],axis=0)
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
        Obsanom = xr.open_dataset('./Observations/ERA5_Meridional_Mean_Anomaly.nc')
        Obsanom = interpolate_obs(Obsanom, DS_CESM_for['lon'])
        OLR_anom = Obsanom['olr'].to_dataset().rename({'olr':yml_usr_info['forecast_olr_name']})
        U200_anom = Obsanom['uwnd200'].to_dataset().rename({'uwnd200':yml_usr_info['forecast_u200_name']})
        U850_anom = Obsanom['uwnd850'].to_dataset().rename({'uwnd850':yml_usr_info['forecast_u850_name']})

        for enen in range(nensembs):
            ### OLR anomaly filtering:
            tmpREolr=OLR_anom.sel(time=slice(first_date_120,first_date))
            tmpREolr=tmpREolr.drop('dayofyear')
            fused_RE_for_OLR = xr.concat([tmpREolr,OLR_cesm_anom.sel(ensemble=enen).to_dataset()],dim='time')
            fused_RE_for_OLR_rolled = fused_RE_for_OLR.rolling(time=120, center=False,min_periods=1).mean().sel(time=slice(OLR_cesm_anom.time.values[0],OLR_cesm_anom.time.values[-1]))
            OLR_cesm_anom_filterd[enen,:,:]=OLR_cesm_anom.sel(ensemble=enen).values - fused_RE_for_OLR_rolled[olrv].values

            ### U200 anomaly filtering:
            tmpRE200=U200_anom.sel(time=slice(first_date_120,first_date))
            tmpRE200=tmpRE200.drop('dayofyear')
            fused_RE_for_200 = xr.concat([tmpRE200,U200_cesm_anom.sel(ensemble=enen).to_dataset()],dim='time')
            fused_RE_for_200_rolled = fused_RE_for_200.rolling(time=120, center=False,min_periods=1).mean().sel(time=slice(U200_cesm_anom.time.values[0],U200_cesm_anom.time.values[-1]))
            U200_cesm_anom_filterd[enen,:,:]=U200_cesm_anom.sel(ensemble=enen).values - fused_RE_for_200_rolled[u200v].values

            ### U850 anomaly filtering:
            tmpRE850=U850_anom.sel(time=slice(first_date_120,first_date))
            tmpRE850=tmpRE850.drop('dayofyear')
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


        MJO_for.attrs["title"] = "CESM2 - MJO RMM Forecast eof(u850,u200,OLR)"
        MJO_for.attrs["description"] = "MJO Forecast in the Prescribed Forecast dataset calculated as in Wheeler and Hendon 2004, a 120-day filter, 15S-15N averaged variables "
        MJO_for.attrs["notes"] = " ONLY Variables RMM1 and RMM2 have been flipped and switched -from eofs_save- to match standard MJO conventions"

        MJO_for.RMM1.attrs['units'] = 'stddev'
        MJO_for.RMM1.attrs['standard_name']='RMM1'
        MJO_for.RMM1.attrs['long_name']='RMM1'

        MJO_for.RMM2.attrs['units'] = 'stddev'
        MJO_for.RMM2.attrs['standard_name']='RMM2'
        MJO_for.RMM2.attrs['long_name']='RMM2'

        MJO_for.OLR_norm.attrs['units'] = 'stddev'
        MJO_for.OLR_norm.attrs['standard_name']='OLR - normalized'
        MJO_for.OLR_norm.attrs['long_name']='OLR - normalized'

        MJO_for.u200_norm.attrs['units'] = 'stddev'
        MJO_for.u200_norm.attrs['standard_name']='u200 normalized'
        MJO_for.u200_norm.attrs['long_name']='u200 normalized'

        MJO_for.u850_norm.attrs['units'] = 'stddev'
        MJO_for.u850_norm.attrs['standard_name']='u200 normalized'
        MJO_for.u850_norm.attrs['long_name']='u200 normalized'

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


    def create_forecasts(self):
        """
        Function to create forecast files for a given range of latitude and settings.

        Parameters:
            yml_data (dict): Data dictionary containing forecast information.
            lons_forecasts (list): List of longitudes for the forecast.

        Returns:
            DS_CESM_for (xr.Dataset): Forecast dataset.
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
        FN_Uwind = sorted(glob.glob(yml_usr_info['forecast_data_loc'] + yml_usr_info['forecast_data_name_str']))
        if len(FN_Uwind) == 0:
            raise FileNotFoundError(f"Files '{yml_usr_info['forecast_data_name_str']}' do not exist... "
                                    f"check your datestring of the filenames.")

        #get the driver dataframe
        DF_Uwind = make_DF_ense(FN_Uwind)

        #get the climatology to create the anomalies
        DS_climo_forecast = self.get_forecast_LT_climo(self.yml_data, self.forecast_lons)

        #loop to make each forecast file:
        for FileCounter, eee in enumerate(range(0, len(DF_Uwind))):
            svname = yml_usr_info['output_files_loc'] + yml_usr_info['output_files_string'] + '_' + DF_Uwind['Init'][eee] + '.nc'

            if os.path.exists(svname):
                print(svname)
                print('The above forecast file already exists... Im skipping it... erase it if you want to make it again')
                continue  # Skip to the next iteration if the file exists

            #Load check and adjust the forecast file and get key variables
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
                donunt=1
                if yml_usr_info['use_forecast_climo']:
                    print('im using the LT dp climo')
                    U850_cesm_anom,U200_cesm_anom,OLR_cesm_anom = self.anomaly_LTD(self.yml_data,DS_CESM_for,DS_climo_forecast,numdays_out)
                else: 
                    U850_cesm_anom,U200_cesm_anom,OLR_cesm_anom = self.anomaly_ERA5(self.yml_data,DS_CESM_for,DS_climo_forecast,numdays_out)
                #function to run the anomaly... get_forecast_anom(yml_data,DS_CESM_for,DS_climo_forecast)
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
            

        return DS_CESM_for,OLR_cesm_anom_filterd,U200_cesm_anom_filterd,U850_cesm_anom_filterd
