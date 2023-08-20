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
import sys
from MJOcast.utils.WHtools import interpolate_obs
from MJOcast.utils.WHtools import check_or_create_paths

class MJOobsProcessor:
    """An example docstring for a class definition."""
    
    def __init__(self,yaml_file_path,scaling_dict=None):
        
        #Global attributes that the functions need! 
        
        with open(yaml_file_path, 'r') as file:
            yml_data = yaml.safe_load(file)
        yml_data
        
        self.yml_data = yml_data
        self.yml_usr_info = yml_data['user_defined_info']
        DSforexample = check_or_create_paths(yml_data)
        self.forecast_lons = DSforexample['lon']
        self.base_dir = yml_data['base_dir']
        self.scaling_dict = scaling_dict
        pass
    
    
    def save_out_obs(self,tot_dict, u200, u850, olr):
        """
        Save the observed MJO dataset to a NetCDF file.

        Parameters:
            tot_dict (dict): Dictionary containing MJO phase, RMM indices, and EOFs.
            u200 (xarray.DataArray): Xarray DataArray containing normalized u200 data.
            u850 (xarray.DataArray): Xarray DataArray containing normalized u850 data.
            olr (xarray.DataArray): Xarray DataArray containing normalized OLR data.

        Returns:
            MJO_fobs (xarray.Dataset): Xarray dataset containing the observed MJO data.
        """

        # Create a new xarray dataset to store the observed MJO data
        MJO_fobs = xr.Dataset(
            {
                "RMM1_obs": (["time"], tot_dict['RMM1_obs']),
                "RMM2_obs": (["time"], tot_dict['RMM2_obs']),
                "RMMind_obs": (["time"], tot_dict['RMMind']),
                "RMMphase_obs": (["time"], tot_dict['MJO_phase']),
                "olr_norm": (["time", "longitude"], olr.data),
                "eof1_olr": (["longitude"], tot_dict['eof1_olr']),
                "eof2_olr": (["longitude"], tot_dict['eof2_olr']),
                "eof1_u850": (["longitude"], tot_dict['eof1_u850']),
                "eof2_u850": (["longitude"], tot_dict['eof2_u850']),
                "eof1_u200": (["longitude"], tot_dict['eof1_u200']),
                "eof2_u200": (["longitude"], tot_dict['eof2_u200']),
                "u200_norm": (["time", "longitude"], u200.data),
                "u850_norm": (["time", "longitude"], u850.data),
            },
            coords={
                "time": olr.time,
                "longitude": olr.lon.values,
            },
        )

        # Add attributes to the dataset
        MJO_fobs.attrs["title"] = "MJO RMM Forecast eof(u850,u200,olr)"
        MJO_fobs.attrs["description"] = "MJO obs in the dataset calculated as in Wheeler and Hendon 2004, a 120-day filter, 15S-15N averaged variables "
        MJO_fobs.attrs["author"] = "S2S_WH_MJO_Forecast_Research_Toolbox"
        MJO_fobs.attrs["questions"] = "wchapman@ucar.edu"

        MJO_fobs.RMM1_obs.attrs['units'] = 'stddev'
        MJO_fobs.RMM1_obs.attrs['standard_name'] = 'RMM1'
        MJO_fobs.RMM1_obs.attrs['long_name'] = 'RMM1'

        MJO_fobs.RMM2_obs.attrs['units'] = 'stddev'
        MJO_fobs.RMM2_obs.attrs['standard_name'] = 'RMM2'
        MJO_fobs.RMM2_obs.attrs['long_name'] = 'RMM2'

        MJO_fobs.olr_norm.attrs['units'] = 'stddev'
        MJO_fobs.olr_norm.attrs['standard_name'] = 'outgoing longwave - normalized'
        MJO_fobs.olr_norm.attrs['long_name'] = 'outgoing longwave - normalized'

        MJO_fobs.u200_norm.attrs['units'] = 'stddev'
        MJO_fobs.u200_norm.attrs['standard_name'] = 'u200 normalized'
        MJO_fobs.u200_norm.attrs['long_name'] = 'u200 normalized'

        MJO_fobs.u850_norm.attrs['units'] = 'stddev'
        MJO_fobs.u850_norm.attrs['standard_name'] = 'u200 normalized'
        MJO_fobs.u850_norm.attrs['long_name'] = 'u200 normalized'

        # Load the BOM index data
        print('...attaching the BOM index for verification...')
        #dataNOAA = pd.read_csv(os.path.join(self.base_dir,'./Observations/BOM_INDEX.txt'), delimiter='\s+', header=0)
        #dataNOAA.columns = ["year", "month", "day", "RMM1", "RMM2", "phase", 'amplitude', 'doggy']

        # Create new variables for the BOM index in the MJO dataset
        #MJO_fobs['RMM1_obs_BOM'] = copy.deepcopy(MJO_fobs['RMM1_obs'].squeeze())
        #MJO_fobs['RMM2_obs_BOM'] = copy.deepcopy(MJO_fobs['RMM2_obs'].squeeze())

        # Loop through the time dimension and assign BOM index values to the MJO dataset
        #for ll in range(len(MJO_fobs['time'])):
        #    yrnum = int(MJO_fobs['time.year'][ll])
        #    monum = int(MJO_fobs['time.month'][ll])
        #    daynum = int(MJO_fobs['time.day'][ll])
        #    ind_date = dataNOAA[(dataNOAA['year'] == yrnum) & (dataNOAA['month'] == monum) & (dataNOAA['day'] == daynum)].index[0]
        #    MJO_fobs['RMM1_obs_BOM'][ll:ll+1] = np.array(dataNOAA[ind_date:ind_date+1]['RMM1'])
        #    MJO_fobs['RMM2_obs_BOM'][ll:ll+1] = np.array(dataNOAA[ind_date:ind_date+1]['RMM2'])

        # Save the MJO dataset to a NetCDF file
        svname = os.path.join(self.base_dir,self.yml_usr_info['obs_data_loc'] + 'MJO_obs_created.nc')
        MJO_fobs.to_netcdf(svname)

        return MJO_fobs
    

    def check_MJO_orientation(self,eof_list, pcs, lons):
        """
        Check the orientation of MJO's first two Empirical Orthogonal Functions (EOFs).

        Parameters:
            eof_list (list): List of empirical orthogonal functions (EOFs).
            pcs (xarray.Dataset): Xarray dataset containing the principal components (PCs).
            lons (array_like): Longitudes.

        Returns:
            loc1 (int): Index of the first dominant EOF.
            loc2 (int): Index of the second dominant EOF.
            scale1 (int): Scaling factor for the first dominant EOF.
            scale2 (int): Scaling factor for the second dominant EOF.
        """
        if self.scaling_dict is not None:
            self.pass_eof_scaling_factors(scaling_dictionary=self.scaling_dict)
            loc1 = self.loc1
            loc2 = self.loc2 
            scale1 = self.scale1
            scale2 = self.scale2 
            return loc1, loc2, scale1, scale2
        else:    

            # Get the first and second EOFs for OLR
            eof1_olr = eof_list[0][0, :]
            eof2_olr = eof_list[0][1, :]
            #print('len lons:',len(lons))
            #print(np.max(np.abs(eof1_olr)),'max abs1!')
            #print(np.max(np.abs(eof2_olr)),'max abs2!')
            # Find the longitude indices of maximum values for the first and second EOFs
            maxolr1_loc = int(np.where(np.abs(eof1_olr.squeeze()) == np.max(np.abs(eof1_olr)))[0])
            maxolr2_loc = int(np.where(np.abs(eof2_olr.squeeze()) == np.max(np.abs(eof2_olr)))[0])

            #print(maxolr1_loc,'loc max abs1!')
            #print(maxolr2_loc,'loc max abs2!')
            #print(np.abs(eof1_olr.squeeze()),'checkthis')

            # Check the orientation of MJO's first two EOFs
            if maxolr1_loc > maxolr2_loc:
                loc1 = 0
                loc2 = 1
            else:
                loc1 = 1
                loc2 = 0

            # Determine the scaling factors for the first two EOFs based on their signs

            if loc1 ==0:
                if eof1_olr[maxolr1_loc] > 0:
                    scale1 = -1
                else:
                    scale1 = 1

                if eof2_olr[maxolr2_loc] > 0:
                    scale2 = 1
                else:
                    scale2 = -1

            elif loc1 ==1:
                if eof1_olr[maxolr1_loc] > 0:
                    scale1 = 1
                else:
                    scale1 = -1

                if eof2_olr[maxolr2_loc] > 0:
                    scale2 = -1
                else:
                    scale2 = 1


            self.loc1 = loc1
            self.loc2 = loc2 
            self.scale1 = scale1
            self.scale2 = scale2 
            return loc1, loc2, scale1, scale2


    def get_phase_and_eofs(self, eof_list, pcs, lons):
        """
        Calculate MJO phase and related EOFs.

        Parameters:
            eof_list (list): List of empirical orthogonal functions (EOFs).
            pcs (xarray.Dataset): Xarray dataset containing the principal components (PCs).
            lons (array_like): Longitudes.

        Returns:
            tot_dict (dict): Dictionary containing MJO phase, RMM indices, and EOFs.
        """

        # Get the locations and scales for MJO EOF1 and EOF2
        loc1, loc2, scale1, scale2 = self.check_MJO_orientation(eof_list, pcs, lons)

        # Calculate scaled EOFs for OLR, U850, and U200 for MJO EOF1 and EOF2
        eof1_olr = eof_list[0][loc1, :] * scale1
        eof2_olr = eof_list[0][loc2, :] * scale2
        eof3_olr = eof_list[0][2, :]

        eof1_u850 = eof_list[1][loc1, :] * scale1
        eof2_u850 = eof_list[1][loc2, :] * scale2
        eof3_u850 = eof_list[1][2, :]

        eof1_u200 = eof_list[2][loc1, :] * scale1
        eof2_u200 = eof_list[2][loc2, :] * scale2
        eof3_u200 = eof_list[2][2, :]

        # Calculate scaled principal components (RMM indices) for MJO EOF1 and EOF2
        pc1 = pcs[:, loc1] * scale1
        pc2 = pcs[:, loc2] * scale2

        # Calculate RMM indices (RMM1_obs and RMM2_obs) based on scaled PCs
        RMM1_obs = pc1
        RMM2_obs = pc2

        # Calculate MJO phase based on RMM indices
        MJO_phase = []
        RMMind = np.sqrt(RMM1_obs**2 + RMM2_obs**2)  # Full index

        for ii in range(RMMind.shape[0]):    
            if np.isnan(RMMind[ii]):
                MJO_phase.append(np.nan)
            elif RMMind[ii] < 1:
                MJO_phase.append(0)
            else:
                ang = np.degrees(np.arctan2(RMM2_obs[ii], RMM1_obs[ii]))
                if ang < 0:
                    ang = ang + 360
                ang = ang + 180

                if ang > 360:
                    ang = ang - 360

                MJO_phase.append(np.floor((ang) / 45) + 1)

        MJO_phase = np.array(MJO_phase)

        # Create a dictionary to store calculated MJO phase, RMM indices, and EOFs
        tot_dict = {
            'MJO_phase': MJO_phase,
            'RMM1_obs': RMM1_obs,
            'RMM2_obs': RMM2_obs,
            'RMMind': RMMind,
            'eof1_olr': eof1_olr,
            'eof2_olr': eof2_olr,
            'eof1_u200': eof1_u200,
            'eof2_u200': eof2_u200,
            'eof1_u850': eof1_u850,
            'eof2_u850': eof2_u850
        }

        return tot_dict

    def make_observed_MJO(self):
        """
        Generate observed MJO indices using observational data.

        Parameters:
            yml_data (dict): Dictionary containing YAML configuration data.
            lons_forecast (array_like): Longitudes of the forecast data.

        Returns:
            OBS_DS (xarray.Dataset): Xarray dataset containing the observed data.
            eof_list (list): List of empirical orthogonal functions (EOFs).
            pcs (xarray.Dataset): Xarray dataset containing the principal components (PCs).
            MJO_fobs (xarray.Dataset): Xarray dataset with observed MJO indices.
        """
        

        # Extract user-defined info from YAML data
        yml_usr_info = self.yml_data['user_defined_info']

        # Open the observed data file based on user settings
        if yml_usr_info['use_era5']:
            obs_fp = os.path.join(self.base_dir,yml_usr_info['obs_data_loc'] + '/ERA5_Meridional_Mean_Anomaly_Filtered120.nc')
            OBS_DS = xr.open_dataset(obs_fp)
        else: 
            print('opening user defined obs file')
            print('-------------------------------------------------')
            print('note: Observed LONS must be from 0->360 and not -180->180')
            print('-------------------------------------------------')
            obs_fp = os.path.join(self.base_dir, yml_usr_info['obs_data_loc'] + '/' + yml_usr_info['usr_named_obs'])
            OBS_DS = xr.open_dataset(obs_fp)

        # Interpolate the observed data to match the forecast longitudes
        OBS_DS = interpolate_obs(OBS_DS, self.forecast_lons)

        # Separate variables from observed dataset for EOF analysis
        OLR_anom2 = OBS_DS['olr'].to_dataset()
        U850_anom2 = OBS_DS['uwnd850'].rename('uwnd').to_dataset()
        U200_anom2 = OBS_DS['uwnd200'].rename('uwnd').to_dataset()

        olr = OLR_anom2['olr']
        u850 = U850_anom2['uwnd']
        u200 = U200_anom2['uwnd']

        # Preprocessing for principal component analysis (PCA)
        olr_var = olr.var("time")
        u850_var = u850.var("time")
        u200_var = u200.var("time")

        olr_norm = np.sqrt(olr_var.mean("lon"))
        u850_norm = np.sqrt(u850_var.mean("lon"))
        u200_norm = np.sqrt(u200_var.mean("lon"))

        olr = olr / olr_norm
        u850 = u850 / u850_norm
        u200 = u200 / u200_norm
        olr[:, -1] = 0
        u200[:, -1] = 0
        u850[:, -1] = 0

        # Perform Empirical Orthogonal Function (EOF) analysis
        print(' ----- Taking the EOF ----- ')
        solver = MultivariateEof([np.array(olr), np.array(u850), np.array(u200)])
        varfrac = solver.varianceFraction()
        
        self.varfrac = varfrac

        # Plot and display variance explained by each EOF
        self.plot_varfrac(varfrac)

        # Extract EOFs and PCs
        eof_list = solver.eofs()
        pcs = solver.pcs(pcscaling=1)

        # Plot observed EOFs and PCs
        self.plot_obs_eof(eof_list, pcs, varfrac, OBS_DS['lon'])

        # Check MJO orientation based on EOFs and PCs
        loc1, loc2, scale1, scale2 = self.check_MJO_orientation(eof_list, pcs, OBS_DS['lon'])

        eof_dict = {'solver':solver,'olr_norm':olr_norm,'u850_norm':u850_norm,'u200_norm':u200_norm,
                   'eof_list':eof_list,'loc1':loc1,'loc2':loc2,'scale1':scale1,'scale2':scale2,'pcs':pcs,
                   'eof_list':eof_list}

        # Get phase and EOFs based on orientation
        tot_dict = self.get_phase_and_eofs(eof_list, pcs, OBS_DS['lon'])
        print('...done making observed EOFS, check ./output_plots/*.png for verification metrics...')

        # Save observed MJO indices to a dataset
        MJO_fobs = self.save_out_obs(tot_dict, u200, u850, olr)
        print('...saved OLR obs file...')
        
        self.OBS_DS = OBS_DS
        self.eof_list = eof_list
        self.pcs = pcs 
        self.MJO_fobs = MJO_fobs
        self.eof_dict = eof_dict
        self.obs_lons = OBS_DS['lon']
        self.varfrac = varfrac

        return OBS_DS, eof_list, pcs, MJO_fobs, eof_dict

    def plot_varfrac(self,varfrac):
        """
        Plot the fraction of the total variance represented by each EOF.

        Parameters:
            varfrac (numpy.array): Array containing the variance fraction values for each EOF.
        """

        # Create a bar plot to show variance fraction for each EOF
        fig = plt.figure(figsize=(15, 10))
        fig.suptitle('Multivariate EOF of OBS (OLR, U850, U200 - 15S-15N average)', fontsize=14)
        eof_num = range(1, 16)
        plt.bar(eof_num, varfrac[0:15], width=0.5)
        plt.axhline(0, color='k')
        plt.xticks(range(1, 16))
        plt.title('Fraction of the total variance represented by each EOF')
        plt.xlabel('EOF #')
        plt.ylabel('Variance Fraction')
        plt.xlim(1, 15)
        plt.ylim(np.min(varfrac), np.max(varfrac) + 0.01)

        # Save the plot as an image
        plt.savefig(self.yml_usr_info['output_plot_loc'] + '/'+'/observed_variance_fraction.png') 
        plt.close()

    def plot_eof(self,ax):
        """
        Customize the EOF plot axis.

        Parameters:
            ax (matplotlib.axes.Axes): The axis of the plot to customize.

        Returns:
            ax (matplotlib.axes.Axes): The customized axis of the plot.
        """

        plt.xlabel('Longitude')
        positions = np.arange(0, 365, 15)
        labels = np.arange(0, 365, 15)
        plt.xticks(positions)
        plt.ylabel('')
        plt.legend()

        return ax

    def plot_obs_eof(self,eof_list, pcs, varfrac, lons):
        """
        Plot the spatial structures of EOF1, EOF2, and EOF3.

        Parameters:
            eof_list (list): List of EOF arrays for OLR, U850, and U200.
            pcs (numpy.array): Array containing principal components.
            varfrac (numpy.array): Array containing the variance fraction values for each EOF.
            lons (numpy.array): Array of longitudes.

        Returns:
            None
        """

        loc1, loc2, scale1, scale2 = self.check_MJO_orientation(eof_list, pcs, lons)

        # Calculate scaled EOF arrays and principal components
        eof1_olr = eof_list[0][loc1, :] * scale1
        eof2_olr = eof_list[0][loc2, :] * scale2
        eof3_olr = eof_list[0][2, :]

        eof1_u850 = eof_list[1][loc1, :] * scale1
        eof2_u850 = eof_list[1][loc2, :] * scale2
        eof3_u850 = eof_list[1][2, :]

        eof1_u200 = eof_list[2][loc1, :] * scale1
        eof2_u200 = eof_list[2][loc2, :] * scale2
        eof3_u200 = eof_list[2][2, :]

        pc1 = pcs[:, loc1] * scale1
        pc2 = pcs[:, loc2] * scale2
        pc3 = pcs[:, 2] * 1

        # Plot EOF1
        fig = plt.figure(figsize=(15, 15))
        fig.suptitle('Spatial structures of EOF1 and EOF2', fontsize=14)

        ax = fig.add_subplot(311)
        plt.title('EOF1 (' + str(int(varfrac[0] * 100)) + '%)', fontsize=10)
        plt.plot(eof1_olr, color='k', linewidth=2, linestyle='solid', label='OLR')
        plt.plot(eof1_u850, color='r', linewidth=2, linestyle='dashed', label='U850')
        plt.plot(eof1_u200, color='b', linewidth=2, linestyle='dotted', label='U200')
        plt.axhline(0, color='k')
        self.plot_eof(ax)

        # Plot EOF2
        ax = fig.add_subplot(312)
        plt.title('EOF2 (' + str(int(varfrac[1] * 100)) + '%)', fontsize=10)
        plt.plot(eof2_olr, color='k', linewidth=2, linestyle='solid', label='OLR')
        plt.plot(eof2_u850, color='r', linewidth=2, linestyle='dashed', label='U850')
        plt.plot(eof2_u200, color='b', linewidth=2, linestyle='dotted', label='U200')
        plt.axhline(0, color='k')
        self.plot_eof(ax)

        # Plot EOF3
        ax = fig.add_subplot(313)
        plt.title('EOF3 (' + str(int(varfrac[2] * 100)) + '%)', fontsize=10)
        plt.plot(eof3_olr, color='k', linewidth=2, linestyle='solid', label='OLR')
        plt.plot(eof3_u850, color='r', linewidth=2, linestyle='dashed', label='U850')
        plt.plot(eof3_u200, color='b', linewidth=2, linestyle='dotted', label='U200')
        plt.axhline(0, color='k')
        self.plot_eof(ax)

        # Save the plot as an image
        plt.savefig(self.yml_usr_info['output_plot_loc'] + '/'+'./observed_eofs.png')
        plt.close()
        
        
    def check_obs_eofs(self):
        """
        Checks the correlation between observed dataset EOF values and ERA5 EOF values for specific modes.

        This function opens an observed dataset file, interpolates the data, and calculates the correlation coefficients
        between specific EOF modes of the observed dataset and corresponding ERA5 EOF values. It compares the calculated
        correlations to a predefined threshold and raises an assertion error if any correlation is below the threshold.

        Returns:
            bool: True if all correlations are above the threshold, indicating a high correlation between observed and ERA5 EOF values.
        """
        fp = self.base_dir + '../../tests/test_cases/eofs_MJO.nc'
        check_corr_ds = xr.open_dataset(fp)
        check_corr_ds = interpolate_obs(check_corr_ds, self.forecast_lons)
        corr_dict = {}

        check_modes = ['eof1_olr', 'eof2_olr', 'eof1_u200', 'eof2_u200', 'eof1_u850', 'eof2_u850']

        for cm in check_modes:
            corrnum = np.corrcoef(check_corr_ds[cm].values, self.MJO_fobs[cm])[0, 1]
            corr_dict[cm] = corrnum

        print('the correlation of the observed dataset EOF values and the ERA5 EOF values is:')
        print(corr_dict)

        for cm in check_modes:
            assert corr_dict[cm] > 0.7, (
            "One of your modes is not correlated highly (<0.7) with the ERA5 observations,\n"
            "consider providing your own scaling factors (ProObs.MJOobsProcessor(yaml_file_path, scaling_dict=scaling_dict)\n"
            "see Wheeler and Hendon 2004 to make sure they are right)\n"
            "If you did pass a scaling dictionary... it is likely not correct"
            )
            
        self.passes_obs_correlation_test = True

        return self.passes_obs_correlation_test
    
    def pass_eof_scaling_factors(self, scaling_dictionary=None):
        """
        Sets the scaling factors for MJO observation EOF modes based on a user-defined dictionary.

        This function allows users to provide a dictionary containing scaling factors for EOF modes.
        If no scaling dictionary is provided, default scaling factors are set and determined for the observations.

        Args:
            scaling_dict (dict, optional): A dictionary containing scaling factors for EOF modes.
                The dictionary should include 'loc1', 'loc2', 'scale1', and 'scale2' keys, representing
                the location of the EOF (0/1) and scaling factors (1 or -1) for the EOF modes 
                1 & 2 in the WH RMM calculation.

        Raises:
            KeyError: If the provided scaling_dict is missing any of the required keys.
        """
        if self.scaling_dict is None:
            self.loc1 = 0
            self.loc2 = 1
            self.scale1 = 1
            self.scale2 = 1
        else:
            try:
                self.loc1 = scaling_dictionary['loc1']
                self.loc2 = scaling_dictionary['loc2']
                self.scale1 = scaling_dictionary['scale1']
                self.scale2 = scaling_dictionary['scale2']
            except KeyError:
                print("The scaling dictionary must contain 'loc1', 'loc2', 'scale1', and 'scale2' keys.")