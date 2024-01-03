import yaml
import os
import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import copy
import re
from datetime import datetime

def func1(arg1,arg2):
    """
    This function takes two arguements setting the first to the second
    
    :param arg1:
    :param arg2:
    :return:
    """
    
    arg1 = arg2
    
def Create_Driver_Yaml(filename, user_, base_dir, use_era5, usr_named_obs, obs_data_loc, forecast_data_loc, forecast_data_name_str, 
                forecast_olr_name, forecast_u200_name, forecast_u850_name, forecast_ensemble_dimension_name, 
                output_plot_loc, output_files_loc, output_files_string, use_forecast_climo, use_observed_climo, 
                regenerate_climo, use_dask_for_climo):
    """
    Create a YAML configuration file from the provided parameters.

    Parameters:
    filename (str): The name of the YAML file to be generated.
    user_ (str): The user name.
    base_dir (str): The base directory.
    use_era5 (bool): Whether to use ERA5 data.
    usr_named_obs (str): Alternative name for observations file.
    obs_data_loc (str): Location of observation data.
    forecast_data_loc (str): Location of forecast data.
    forecast_data_name_str (str): String pattern for forecast data filenames.
    forecast_olr_name (str): Name of OLR forecast variable.
    forecast_u200_name (str): Name of 200mb wind forecast variable.
    forecast_u850_name (str): Name of 850mb wind forecast variable.
    forecast_ensemble_dimension_name (str): Name of the ensemble dimension.
    output_plot_loc (str): Location for output plots.
    output_files_loc (str): Location for output files.
    output_files_string (str): String for naming MJO forecast files.
    use_forecast_climo (bool): Whether to use forecast climatology.
    use_observed_climo (bool): Whether to use observed climatology.
    regenerate_climo (bool): Whether to regenerate climatology.
    use_dask_for_climo (bool): Whether to use Dask for climatology processing.
    """
    Total_yaml = {
        'user': user_,
        'base_dir': base_dir,
        'user_defined_info': {
            'use_era5': use_era5,
            'usr_named_obs': usr_named_obs,
            'obs_data_loc': obs_data_loc,
            'forecast_data_loc': forecast_data_loc,
            'forecast_data_name_str': forecast_data_name_str,
            'forecast_olr_name': forecast_olr_name,
            'forecast_u200_name': forecast_u200_name,
            'forecast_u850_name': forecast_u850_name,
            'forecast_ensemble_dimension_name': forecast_ensemble_dimension_name,
            'output_plot_loc': output_plot_loc,
            'output_files_loc': output_files_loc,
            'output_files_string': output_files_string,
            'use_forecast_climo': use_forecast_climo,
            'use_observed_climo': use_observed_climo,
            'regenerate_climo': regenerate_climo,
            'use_dask_for_climo': use_dask_for_climo,
            'regenerate_climo': regenerate_climo
        }
    }

    with open(filename, 'w') as yaml_file:
        yaml.dump(Total_yaml, yaml_file, default_flow_style=False)

# Call the function with the desired variable values and filename
#Create_Yaml('output.yaml', user_, base_dir, use_era5, usr_named_obs, obs_data_loc, forecast_data_loc, forecast_data_name_str, 
#            forecast_olr_name, forecast_u200_name, forecast_u850_name, forecast_ensemble_dimension_name, 
#            output_plot_loc, output_files_loc, output_files_string, use_forecast_climo, use_observed_climo, 
#            regenerate_climo, use_dask_for_climo)


def flip_lat_if_necessary(data):
    """
    Check the orientation of the latitude dimension in an xarray dataset and flip it if necessary.

    Parameters:
        data (xr.Dataset): Input xarray dataset.

    Returns:
        data_flipped (xr.Dataset): Flipped xarray dataset, if necessary.
    """
    # Get the latitude values
    lat_values = data.coords["lat"].values
    
    # Check if latitude is oriented from South to North
    if lat_values[0] < lat_values[-1]:
        # Latitude is oriented from South to North, no need to flip
        return data
    else:
        # Latitude is oriented from North to South, flip it
        data_flipped = data.reindex(lat=lat_values[::-1])
        return data_flipped

def switch_lon_to_0_360(data):
    """
    Check if the longitude values in the xarray dataset are in the range -180 to 180 degrees,
    and switch them to the range 0 to 360 degrees if needed.

    Parameters:
        data (xr.Dataset or xr.DataArray): Input xarray dataset or data array.

    Returns:
        data_with_lon_0_360 (xr.Dataset or xr.DataArray): Xarray dataset or data array with longitude values
        in the range 0 to 360 degrees.
    """
    # Check if the data is a DataArray
    is_data_array = isinstance(data, xr.DataArray)

    if is_data_array:
        lon_var = data.coords['lon']
        if lon_var.min() >= 0 and lon_var.max() <= 360:
            # The data already has longitude values in the range 0 to 360 degrees
            return data
        else:
            # Switch the longitude values to the range 0 to 360 degrees
            data.coords['lon'] = (lon_var + 360) % 360
            data = data.sortby(data.lon)
            return data
    else:
        # If the data is a Dataset, check for 'lon' or 'longitude' coordinate variables
        lon_var_name = None
        for var_name in data.coords:
            if 'lon' in var_name.lower() or 'longitude' in var_name.lower():
                lon_var_name = var_name
                break

        if lon_var_name is not None:
            lon_var = data.coords[lon_var_name]
            if lon_var.min() >= 0 and lon_var.max() <= 360:
                # The data already has longitude values in the range 0 to 360 degrees
                return data
            else:
                # Switch the longitude values to the range 0 to 360 degrees
                data[lon_var_name] = (lon_var + 360) % 360
                data = data.sortby(data[lon_var_name])
                return data
        else:
            # If 'lon' or 'longitude' coordinate variables are not found, raise an error
            raise ValueError("No 'lon' or 'longitude' coordinate variables found in the dataset.")


def check_lat_lon_coords(data):
    """
    Check if either "latitude/longitude" or "lat/lon" are in the coordinates of the xarray dataset.

    Parameters:
        data (xr.Dataset): Input xarray dataset.

    Returns:
        has_lat_lon_coords (bool): True if either "latitude/longitude" or "lat/lon" are in the coordinates, False otherwise.
    """
    has_lat_lon_coords = False

    # Check for "latitude/longitude" or "lat/lon" in the coordinates
    if "latitude" in data.coords and "longitude" in data.coords:
        has_lat_lon_coords = True
    elif "lat" in data.coords and "lon" in data.coords:
        has_lat_lon_coords = True

    return has_lat_lon_coords


def check_or_create_paths(yml_data):
    """
    Check the paths and files required for the forecast data.

    Parameters:
        yml_data (dict): Dictionary containing information from the YAML file.

    Returns:
        DS (xarray.Dataset): Xarray dataset of the forecast data.
    """

    yml_usr_info = yml_data['user_defined_info']
    check_forecast_data = yml_usr_info['forecast_data_loc']

    try:
        # Check if the forecast data directory exists
        assert os.path.exists(check_forecast_data), print(f"The path '{check_forecast_data}' exists.")

        # Get a list of forecast files using a specified filename pattern
        for_fil_list = sorted(glob.glob(yml_usr_info['forecast_data_loc'] + '/'+ yml_usr_info['forecast_data_name_str']))
        file_count = len(for_fil_list)
        print(f"Number of forecast files to process: {file_count}")

        try:
            # Check if there are any forecast files to process
            assert file_count > 0
        except: 
            raise FileNotFoundError(f"Files '{yml_usr_info['forecast_data_loc'] + '/'+ yml_usr_info['forecast_data_name_str']}' do not exist.")
        # Check the first forecast file and get the xarray dataset
        Bingo, DS = check_forecast_files(for_fil_list[0], yml_usr_info)

        if Bingo:
            print('Initial look at forecast files passes the first test')
        else:
            raise RuntimeError("Did not pass file test, check the printed dictionary to ensure your variable names are defined")

    except AssertionError as e:
        # If the path doesn't exist, raise a custom error message
        raise FileNotFoundError(f"The path '{yml_usr_info['forecast_data_loc'] + '/' +yml_usr_info['forecast_data_name_str']}' does not exist.")

    return DS

def check_forecast_files(for_file_list, yml_usr_info):
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
        raise RuntimeError("it broke while re-orienting your xarray files to be S->N and 0-360 degrees lon") 
    
    # Get variable names from user settings
    u200v = yml_usr_info['forecast_u200_name']
    u850v = yml_usr_info['forecast_u850_name']
    olrv = yml_usr_info['forecast_olr_name']

    # Check if required variables are present in the dataset
    if u200v in DS.variables and u850v in DS.variables and olrv in DS.variables:
        Bingo = True
    else:
        print('it looks like the defined u200,u850,or olr variables are not present in the file')
        Bingo = False

    # Get the ensemble dimension name from user settings
    ensemble_name = yml_usr_info['forecast_ensemble_dimension_name']

    # Check if the ensemble dimension is already present in the dataset
    if ensemble_name in DS.coords:
        print('ensemble dimension length:', len(DS[ensemble_name]))
    else:
        # If the ensemble dimension is not present, try to add it
        try:
            # Create an 'ensemble' coordinate array and assign it to the dataset
            DS = DS.assign_coords({ensemble_name: np.arange(len(DS[ensemble_name]))})
            print('expanding coords to include ensemble')
            print('ensemble dimension length:', len(DS[ensemble_name]))
        except:
            # If there is an error, add 'ensemble' as a new dimension to the dataset
            ensemble_values = [0]
            DS = DS.expand_dims(ensemble=ensemble_values)
            # Assign the coordinate values for the new 'ensemble' dimension
            DS = DS.assign_coords(ensemble=ensemble_values)
            
    if 'time' in DS.coords: 
        print(f"there are {len(DS['time.dayofyear'])} forecast lead days in these files")
    else: 
        raise RuntimeError("the files MUST be compatible with xarray's DS['time.dayofyear'] functionality")

    return Bingo, DS


def convert_dates_to_string(input_string):
    # Define the regex pattern to match various date formats
    date_pattern = r'(\d{1,2}[A-Za-z]{3}\d{4}|\d{1,2}-[A-Za-z]{3}-\d{4}|[A-Za-z]{3}[-_\s]\d{1,2}[-_\s]\d{4}|[0-3]?\d[0-1]?\d\d{2}|[0-1]?\d[0-3]?\d\d{2})'
    
    # Find all matches of the date pattern in the input string
    try:
        matches = re.findall(date_pattern, input_string)
    except:
        raise RuntimeError("the forecast files do not have a proper datestring as defined in the settings.yaml please change their names")
    
    
    # Loop through each matched date and convert it to the desired format
    for match in matches:
        # Parse the date string to a datetime object
        try:
            date_obj = datetime.strptime(match, '%b-%d-%Y')
        except ValueError:
            try:
                date_obj = datetime.strptime(match, '%d-%b-%Y')
            except ValueError:
                try:
                    date_obj = datetime.strptime(match, '%b_%d_%Y')
                except ValueError:
                    try:
                        date_obj = datetime.strptime(match, '%b-%d-%Y')
                    except ValueError:
                        try:
                            date_obj = datetime.strptime(match, '%m%d%y')
                        except ValueError:
                            try:
                                date_obj = datetime.strptime(match, '%d%b%Y')
                            except ValueError:
                                continue  # Skip if the date format doesn't match any of the known formats
        # Convert the datetime object to the desired format
        formatted_date = date_obj.strftime('%d%b%Y')
        
        # Replace the matched date in the input string with the formatted date
        input_string = input_string.replace(match, formatted_date)
    
    return input_string,matches

def interpolate_obs(OBS_DS, lons_forecast):
    """
    Interpolate observed data to match forecast longitudes.

    Parameters:
        OBS_DS (xarray.Dataset): Xarray dataset containing the observed data.
        lons_forecast (array_like): Longitudes of the forecast data.

    Returns:
        OBS_DS (xarray.Dataset): Interpolated xarray dataset of the observed data.
    """

    # Check if the number of forecast longitudes is already 360
    if len(lons_forecast) == 360:
        return OBS_DS
    else:
        for var_name in OBS_DS.coords:
            if 'lon' in var_name.lower() or 'longitude' in var_name.lower():
                lon_var_name = var_name
                break
                
        if lon_var_name=='lon':
            # Interpolate observed data to match forecast longitudes
            OBS_DS = OBS_DS.interp(lon=np.array(lons_forecast)).fillna(0)
        if lon_var_name=='longitude':
            # Interpolate observed data to match forecast longitudes
            OBS_DS = OBS_DS.interp(longitude=np.array(lons_forecast)).fillna(0)

    return OBS_DS

def plot_phase_space(ax):
    """
    Plots a phase space diagram with annotations and labels.
    
    Parameters:
    ax (matplotlib.axes._axes.Axes): The matplotlib axes on which to plot the diagram.
    
    Returns:
    matplotlib.axes._axes.Axes: The axes with the plotted phase space diagram.
    """
    # Set the x and y labels and limits
    plt.xlabel('RMM1', fontsize=14)
    plt.ylabel('RMM2', fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)  # Increase tick label size
    plt.xlim(-4, 4)
    plt.ylim(-4, 4)
    
    # Plot diagonal and cross lines
    plt.plot([-4, 4], [-4, 4], color='k')
    plt.plot([-4, 4], [4, -4], color='k')
    plt.plot([-4, 4], [0, 0], color='k')
    plt.plot([0, 0], [-4, 4], color='k')
    
    # Add a circle to the diagram
    circle1 = plt.Circle((0, 0), 1, color='k', fill=False)
    ax.add_artist(circle1)
    
    # Add text annotations for different regions
    plt.text(-1.4, -3, '2', fontsize=30, color='grey')
    plt.text(1.4, -3, '3', fontsize=30, color='grey')
    plt.text(-0.9, -3.8, 'Indian Ocean', fontsize=20, color='grey')
    plt.text(2.6, -1, '4', fontsize=30, color='grey')
    plt.text(2.6, 1, '5', fontsize=30, color='grey')
    plt.text(3.5, -1, 'Maritime Continent', fontsize=20, color='grey', rotation='vertical')
    plt.text(1.4, 2.6, '6', fontsize=30, color='grey')
    plt.text(-1.4, 2.6, '7', fontsize=30, color='grey')
    plt.text(-0.9, 3.6, 'Pacific Ocean', fontsize=20, color='grey')
    plt.text(-2.6, 1, '8', fontsize=30, color='grey')
    plt.text(-2.6, -1, '1', fontsize=30, color='grey')
    plt.text(-3.8, -1, 'West. Hem. & Africa', fontsize=20, color='grey', rotation='vertical')
    
    # Return the axes with the plotted phase space diagram
    return ax
