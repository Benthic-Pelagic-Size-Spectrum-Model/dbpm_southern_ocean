#!/usr/bin/env python3

# Library of useful functions
# Authors: Denisse Fierro Arcos
# Date last update: 2026-11-04

# Loading libraries
import xarray as xr
import re
import numpy as np
import os
from glob import glob
import pandas as pd
import json
import cftime
from datetime import datetime as dt


#Transforming netCDF files to zarr
def netcdf_to_zarr(file_path, path_out):
    '''
    Inputs:
    - file_path (character) File path where GFDL output is located
    - path_out (character) File path where outputs should be stored as zarr files

    Outputs:
    - None. This function saves results as zarr files in the path provided.
    '''

    #Loading and rechunking data
    da = xr.open_dataarray(file_path).chunk({'lat': '50MB', 'lon': '50MB'})

    #Change date format
    try:
        da['time'] = pd.DatetimeIndex(da.indexes['time'].to_datetimeindex(),
                                      dtype = 'datetime64[ns]')
    except:
        pass

    #Save results
    da.to_zarr(path_out, consolidated = True, mode = 'w')


## Extracting GFDL outputs for region of interest using boolean mask
def extract_gfdl(file_path_or_data_array, mask, path_out, cross_dateline = False):
    '''
    Inputs:
    - file_path_or_data_array (character or Data Array) Provide either file path where GFDL 
    zarr file is located or Data Array
    - mask (boolean data array) Grid cells within region of interest should be identified
    as 1.
    - path_out (character) File path where outputs should be stored as zarr files
    - cross_dateline (boolean) Default is False. If set to True, it will convert longitudes 
    from +/-180 to 0-360 degrees before extracting data for region of interest

    Outputs:
    - None. This function saves results as zarr files in the path provided.
    '''

    #Loading and rechunking data
    if isinstance(file_path_or_data_array, str):
        da = xr.open_zarr(file_path_or_data_array)
        #Getting name of variable contained in dataset
        [var] = list(da.keys())
        da = da[var]
    elif isinstance(file_path_or_data_array, xr.DataArray):
        da = file_path_or_data_array
    else:
        print('Please provide either a file path or data array to run function.')

    #Fix time format if needed
    try:
        new_time = da.indexes['time'].to_datetimeindex()
        da['time'] = new_time
    except:
        pass
    
    #Apply mask and remove rows where all grid cells are empty to reduce data array size
    if cross_dateline:
        da = da.where(mask == 1)
        da['lon'] = da.lon%360
        da = da.sortby('lon')
        da = da.dropna(dim = 'lon', how = 'all').dropna(dim = 'lat', how = 'all')
    else:
        da = da.where(mask == 1, drop = True)
    
    #Rechunking data
    if 'time' in da.dims:
        da = da.chunk({'time': '50MB', 'lat': len(da.lat), 'lon': len(da.lon)})
    else:
        da = da.chunk({'lat': len(da.lat), 'lon': len(da.lon)})

    #Save results
    da.to_zarr(path_out, consolidated = True, mode = 'w')


## Calculating area weighted means
def weighted_mean_timestep(file_paths, weights, region):
    '''
    Inputs:
    - file_paths (list) File paths pointing to zarr files from which weighted means 
    will be calculated and stored in a single data frame
    - weights (data array) Contains the weights to be used when calculating weighted mean.
    It should NOT include NaN, zeroes (0) should be used instead.
    - region (character) Name of the region to be recorded in data frame

    Outputs:
    - df (pandas data frame) containing area weighted mean per time step
    '''
    
    #Loading all zarr files into a single dataset
    da = xr.open_mfdataset(file_paths, engine = 'zarr')
    #Fix time format if needed
    try:
        new_time = da.indexes['time'].to_datetimeindex()
        da['time'] = new_time
    except:
        pass

    #Apply weights
    da_weighted = da.weighted(weights)
    #Calculate weighted mean
    da_weighted_mean = da_weighted.mean(('lat', 'lon'))
    
    #Transform to data frame
    df = da_weighted_mean.to_dataframe().reset_index()

    #Getting names of variables in data frame
    col_names = [i for i in df.columns if i != 'time']
    #Getting name of experiment from file path
    [exp] = re.findall('cobalt2_(.*?)_', file_paths[0])

    #If no depth file is included, add variable and leave it empty
    if 'deptho' not in col_names:
        df['depth_m'] = np.nan
    else:
        col_names.remove('deptho')
        df = df.rename(columns = {'deptho': 'depth_m'})
    
    #Add metadata to data frame
    df['tot_area_m2'] = weights.values.sum()
    df['year'] = df.apply(lambda x: x.time.year, axis = 1)
    df['month'] = df.apply(lambda x: x.time.strftime('%B'), axis = 1)
    df['region'] = region
    df['scenario'] = exp

    #Rearrange columns 
    names = ['region', 'scenario', 'time', 'year', 'month', 
             'depth_m', 'tot_area_m2'] + col_names
    df = df[names]

    return df


#Calculating export ratio
def getExportRatio(folder_gridded_data, gfdl_exp):
    '''
    Inputs:
    - folder_gridded_data (character) File path pointing to folder containing
    zarr files with GFDL data for the region of interest
    - gfdl_exp (character) Select 'ctrl_clim' or 'obs_clim'

    Outputs:
    - sphy (data array) Contains small phytoplankton. This is data frame simply
    renames 'phypico-vint' data
    - lphy (data array) Contains large phytoplankton: difference between 'phyc_vint'
    and 'phypico_vint'
    - er (data array) Contains export ratio
    '''

    #Get list of files in experiment
    file_list = glob(os.path.join(folder_gridded_data, f'*_{gfdl_exp}_*'))
    
    #Load sea surface temperature
    tos = xr.open_zarr([f for f in file_list if '_tos_' in f][0])['tos']
    
    #load depth
    depth = xr.open_zarr([f for f in file_list if '_deptho_' in f][0])['deptho']
    
    #Load phypico-vint
    sphy = xr.open_zarr([f for f in file_list if '_phypico-vint_' in f][0])['phypico-vint']
    #Rename phypico-vint to sphy
    sphy.name = 'sphy'

    #Load phyc-vint
    ptotal = xr.open_zarr([f for f in file_list if '_phyc-vint_' in f][0])['phyc-vint']

    #Calculate large phytoplankton
    lphy = ptotal-sphy
    #Give correct name to large phytoplankton dataset
    lphy.name = 'lphy'

    #Calculate phytoplankton size ratios
    plarge = lphy/ptotal
    psmall = sphy/ptotal

    #Calculate export ration
    er = (np.exp(-0.032*tos)*((0.14*psmall)+(0.74*(plarge)))+
          (0.0228*(plarge)*(depth*0.004)))/(1+(depth*0.004))
    #If values are negative, assign a value of 0
    er = xr.where(er < 0, 0, er)
    #If values are above 1, assign a value of 1
    er = xr.where(er > 1, 1, er)
    er.name = 'export_ratio'
    
    return sphy, lphy, er


#Calculate slope and intercept
def GetPPIntSlope(file_paths, mmin = 10**(-14.25), mmid = 10**(-10.184), 
                  mmax = 10**(-5.25), output = 'both'):
    '''
    Inputs:
    - file_paths (list) File paths pointing to zarr files from which slope and
    intercept will be calculated
    - mmin (numeric)  Default is 10**(-14.25). ????
    - mmid (numeric)  Default is 10**(-10.184). ????
    - mmax (numeric)  Default is 10**(-5.25). ????
    - output (character) Default is 'both'. Select what outputs should be returned. 
    Choose from 'both', 'slope', or 'intercept'

    Outputs:
    - (Data array) - Depends on value of 'output' parameter. 
    '''
    
    #load large phytoplankton
    lphy = xr.open_zarr([f for f in file_paths if '_lphy_' in f][0])['lphy']
    
    #Load small phytoplankton
    sphy = xr.open_zarr([f for f in file_paths if '_sphy_' in f][0])['sphy']
    
    #Convert sphy and lphy from mol C / m^3 to g C / m^3
    sphy = sphy*12.0107
    lphy = lphy*12.0107

    #From Appendix of Barnes 2010 JPR, used in Woodworth-Jefcoats et al 2013
    #GCB. The scaling of biomass with body mass can be described as B=aM^b
    #the exponent b (also equivalent to the slope b in a log B vs log M 
    #relationship) can be assumed:
    #0.25 (results in N ~ M^-3/4) or 0 (results in N ~ M^-1)
    # most studies seem to suggest N~M^-1, so can assume that and test 
    #sensitivity of our results to this assumption. 
      
    #Calculate a and b in log B (log10 abundance) vs. log M (log10 gww)
    #in log10 gww
    midsmall = np.log10((mmin+mmid)/2) 
    midlarge = np.log10((mmid+mmax)/2)

    #convert to log10 (gww/size class median size) for log10 abundance
    small = np.log10((sphy*10)/(10**midsmall))
    #convert to log10 (gww/size class median size) for log10 abundance
    large = np.log10((lphy*10)/(10**midlarge))

    #Calculating lope
    b = (small-large)/(midsmall-midlarge)
    b.name = 'slope'

    #a is really log10(a), same a when small, midsmall are used
    a = large-(b*midlarge)
    a.name = 'intercept'

    # a could be used directly to replace 10^pp in sizemodel()
    if output == 'slope':
        return b
    if output == 'intercept':
        return a
    if output == 'both':
        return a, b


# Calculate spinup from gridded data
def gridded_spinup(file_path, start_spin, end_spin, spinup_period, 
                   mean_spinup = False, **kwargs):
    '''
    Inputs:
    - file_path (character) File path pointing to zarr file from which spinup
    will be calculated
    - start_spin (character or numeric) Year or date spinup starts
    - end_spin (character or numeric) Year or date spinup ends
    - spinup_period (pandas Datetime array) New time labels for spinup period.
    Must be a multiple of spinup range (i.e., difference between start and end
    spin)
    - mean_spinup (boolean) Default is False. If set to True, then the spinup
    period will be based on the mean value over the spin period.
    **Optional**: 
    - file_out (character) File path to save results

    Outputs:
    spinup_da (data array) Spinup data containing information within spinup
    range
    '''
    
    #Loading data
    
    da = xr.open_zarr(file_path)
    #Getting name of variable contained in dataset
    [var] = list(da.keys())
    #Select period to be used for spinup
    da = da[var].sel(time = slice(str(start_spin), str(end_spin)))
    #If spinup should be created based on mean values over spinup period
    if mean_spinup:
        da = da.mean('time')
        spinup_da = [da] * len(spinup_period)
    else:
        spinup_da = [da] * int(len(spinup_period)/len(da.time))

    #Create spinup data array
    spinup_da = xr.concat(spinup_da, dim = 'time')
    spinup_da['time'] = spinup_period
    
    #Updating chunks
    spinup_da = spinup_da.chunk({'time': '50MB', 'lat': 100, 'lon': 240})
    spinup_da = spinup_da.drop_encoding()

    #Save result if path is provided
    if kwargs.get('file_out', None):
        f_out = kwargs.get('file_out')
        #Ensure folder in file path exists
        os.makedirs(os.path.dirname(f_out), exist_ok = True)
        spinup_da.to_zarr(f_out, consolidated = True, mode = 'w')
    
    return spinup_da


#Format gridded_params in Python friendly way
def gridded_param_python(gridded_params):
    '''
    Inputs:
    - gridded_params (dictionary) Containing the output from the `sizeparam` function in R.
    Stored as a `json` file.

    Outputs:
    griddded_python (dictionary) Containing gridded_params in a Python friendly format. New
    entries calculated as needed, and entries not used in the DBPM gridded model were removed.
    '''

    gridded_python = {'timesteps_years': gridded_params['timesteps_years'][0],
                      'numb_time_steps': gridded_params['numb_time_steps'][0],
                      'effort': gridded_params['effort'],
                      'fish_mort_pred': gridded_params['fish_mort_pred'][0],
                      'fish_mort_detritivore': gridded_params['fish_mort_detritivore'][0],
                      'hr_volume_search': gridded_params['hr_volume_search'][0],
                      'detritus_coupling': gridded_params['detritus_coupling'][0],
                      'log10_pred_prey_ratio': gridded_params['log10_pred_prey_ratio'][0],
                      'log_prey_pref': gridded_params['log_prey_pref'][0],
                      'hr_vol_filter_benthos': gridded_params['hr_vol_filter_benthos'][0],
                      'metabolic_req_pred': gridded_params['metabolic_req_pred'][0],
                      'metabolic_req_detritivore': gridded_params['metabolic_req_detritivore'][0],
                      'defecate_prop': gridded_params['defecate_prop'][0],
                      'growth_prop': 1-(gridded_params['defecate_prop'][0]),
                      'def_low': gridded_params['def_low'][0],
                      'high_prop': 1-(gridded_params['def_low'][0]),
                      'growth_pred': gridded_params['growth_pred'][0],
                      'growth_detritivore': gridded_params['growth_detritivore'][0],
                      'growth_detritus': gridded_params['growth_detritus'][0],
                      'energy_pred': gridded_params['energy_pred'][0],
                      'energy_detritivore': gridded_params['energy_detritivore'][0],
                      'handling': gridded_params['handling'][0],
                      'dynamic_reproduction': gridded_params['dynamic_reproduction'][0],
                      'c1': gridded_params['c1'][0],
                      'activation_energy': gridded_params['activation_energy'][0],
                      'boltzmann': gridded_params['boltzmann'][0],
                      'natural_mort': gridded_params['natural_mort'][0],
                      'size_senescence': gridded_params['size_senescence'][0],
                      'exp_senescence_mort': gridded_params['exp_senescence_mort'][0],
                      'const_senescence_mort': gridded_params['const_senescence_mort'][0],
                      'log_size_increase': gridded_params['log_size_increase'][0],
                      'log10_size_bins': gridded_params['log10_size_bins'],
                      'numb_size_bins':  gridded_params['numb_size_bins'][0],
                      'ind_min_pred_size': (gridded_params['ind_min_pred_size'][0])-1,
                      'ind_min_detritivore_size': (gridded_params['ind_min_detritivore_size'][0])-1,
                      'idx_new': list(range(gridded_params['ind_min_detritivore_size'][0],
                                            gridded_params['numb_size_bins'][0])),
                      'ind_min_fish_pred': int(gridded_params['ind_min_fish_pred'][0]-1),
                      'ind_min_fish_det': int(gridded_params['ind_min_fish_det'][0]-1),
                      'idx': (np.array(gridded_params['idx'])-1).tolist(),
                      'init_pred': gridded_params['init_pred'],
                      'init_detritivores': gridded_params['init_detritivores'],
                      'init_detritus': gridded_params['init_detritus'][0]
                     }
    
    return gridded_python


### DPBM functions ----
# Build a lookup table for diet preference. Looks at all combinations of predator 
# and prey body size: diet preference (in the predator spectrum only)
def phi_f(q, log10_pred_prey_ratio, log_prey_pref):
    phi = np.where(q > 0, 
                   np.exp(-(q-log10_pred_prey_ratio)*(q-log10_pred_prey_ratio)/
                          (2*log_prey_pref**2))/(log_prey_pref*np.sqrt(2.0*np.pi)),
                   0) 
    return phi


# Function to build lookup tables for (constant) growth
# Considers components which remain constant
def gphi_f(pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref):
    gphi = 10**(-pred_prey_matrix)*phi_f(pred_prey_matrix, log10_pred_prey_ratio, 
                                         log_prey_pref)
    return gphi


# Function to build lookup tables for (constant) mortality
def mphi_f(rev_pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref, 
           metabolic_req_pred):
    mphi = (10**(metabolic_req_pred*rev_pred_prey_matrix)*
            phi_f(rev_pred_prey_matrix, log10_pred_prey_ratio, log_prey_pref))
    return mphi


# Function to build lookup table for components of 10^(alpha*x)
def expax_f(log10_size_bins, metabolic_req_pred):
    expax = 10**(log10_size_bins*metabolic_req_pred)
    return expax


# Create predator-prey combination matrix
# The square matrix holds the log(predatorsize/preysize) for all combinations of
# predator-prey sizes
def pred_prey_matrix(log10_size_bins):
    '''
    Inputs:
    - log10_size_bins (numpy array) Containing the log of all sizes of predator and
    prey items
    
    Outputs:
    ppm (numpy array) Contains all combinations of predator-prey sizes
    '''
    
    ppm = np.empty([len(log10_size_bins), len(log10_size_bins)])
    for i, j in enumerate(log10_size_bins):
        ppm[:, i] = log10_size_bins[i] - log10_size_bins
    return ppm


# Initialising matrices using size classes and time as dimensions
def init_da(log10_size_bins, time):
    '''
    Inputs:
    - log10_size_bins (numpy array) Containing the log of all sizes of predator and
    prey items
    - time (numpy array) Containing dates to be included in final data array
    
    Outputs:
    da (data array) Containing zeroes with size classes and time as dimensions
    '''
    data_start = np.zeros((len(log10_size_bins), len(time)))
    da = xr.DataArray(data = data_start, dims = ['size_class', 'time'], 
                      coords = {'size_class': log10_size_bins, 'time': time})
    return da


# Gravity model ----
# Redistribute total effort across grid cells according to proportion of biomass in that 
# grid cell using graivity model, Walters & Bonfil, 1999, Gelchu & Pauly 2007 ideal free
# distribution - Blanchard et al 2008
def gravitymodel(effort, prop_b, depth, n_iter, sea_ice_mask):
    '''
    Inputs:
    - effort (Data array) Fishing effort for a single time step
    - prop_b (Data array) Proportion of total fishable biomass for each grid cell at a 
    single time step
    - depth (Data array) Bathymetry of the area of interest
    - n_iter (integer) Number of iterations needed to redistribute fishing effort
    - sea_ice_mask (Data array) Mask describing sea ice covered areas inaccessible to 
    fishers
    
    Outputs:
    eff (Data array) Containing redistributed fishing effort
    '''

    eff = effort
    d = (1-(depth/depth.max()))

    #Initialise loop
    i = 1
    while(i <= n_iter):
        suit = prop_b*d*sea_ice_mask
        rel_suit = suit/(suit.sum())
        neweffort = eff+(rel_suit*eff)
        mult = (eff.sum())/(neweffort.sum())
        eff = neweffort*mult
        i += i  

    #Ensuring NA values are transformed to zeroes
    eff = xr.where(np.isfinite(eff), eff, 0)
    return eff


# Calculate fishing mortality and effort ------
def effort_calculation(predators, detritivores, effort, depth, size_bin_vals, 
                       gridded_params, sea_ice_mask):
    '''
    Inputs:
    - predators (2D data array) Pelagic predator biomass
    - detritivores (2D data array) Bethic detritivore biomass
    - effort (2D data array) Fishing effort
    - depth (2D data array) Bathymetry of the area of interest
    - size_bins_vals (1D data array) Size classes in grams
    - gridded_params (dictionary) DBPM parameters
    - sea_ice_mask (Data array) Mask describing sea ice covered areas inaccessible
    to fishers

    Outputs:
    - new_effort (2D data array) Fishing effort calculated for next time step
    '''
    if effort.sum() != 0:
        pred_bio = ((predators*size_bin_vals*gridded_params['log_size_increase']).
                    isel(size_class = slice(gridded_params['ind_min_fish_pred'],
                                            -1))).sum('size_class')
    
        det_bio = ((detritivores*size_bin_vals*gridded_params['log_size_increase']).
                   isel(size_class = slice(gridded_params['ind_min_fish_det'],
                                           -1))).sum('size_class')
    
        sum_bio = pred_bio+det_bio
    
        prop_b = sum_bio/sum_bio.sum()
    
        #Calculate new effort
        new_effort = gravitymodel(effort, prop_b, depth, 1, sea_ice_mask)
        new_effort = xr.where(np.isfinite(depth), new_effort, np.nan)
    else:
        new_effort = effort
    
    # Adjusting time stamp
    try:
        new_effort['time'] = effort.time.values
    except:
        new_effort['time'] = [effort.time.values]
    try:
        new_effort = new_effort.transpose('time', 'lat', 'lon')
    except:
        pass
    #Adding name
    new_effort.name = 'effort'

    #Return effort
    return new_effort


# Loading fixed gridded DBPM inputs ------
def loading_dbpm_fixed_inputs(base_folder):
    '''
    Inputs:
    - base_folder (character) Full path to folder where fixed DBPM gridded inputs
    are stored

    Outputs:
    - dbpm_inputs (xarray Dataset) Contains fixed gridded inputs needed to run DBPM
    '''
    #Loading data from base folder
    #Preference benthic
    pref_benthos = xr.open_zarr(glob(
        os.path.join(base_folder, 
                     'pref-benthos_*'))[0])['pref_benthos']
    #Preference pelagic
    pref_pelagic = xr.open_zarr(glob(
        os.path.join(base_folder, 
                     'pref-pelagic_*'))[0])['pref_pelagic']
    #Constant growth
    constant_growth = xr.open_zarr(glob(
        os.path.join(base_folder, 
                     'const-growth_*'))[0])['constant_growth']
    #Constant mortality
    constant_mortality = xr.open_zarr(glob(
        os.path.join(base_folder, 
                     'const-mort_*'))[0])['constant_mortality']
    #Consumption pelagic
    consume_pelagic = xr.open_zarr(glob(
        os.path.join(base_folder, 
                     'consume-pelagic_*'))[0])['consume_pelagic']
    #Consumption benthos
    consume_benthos = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'consume-benthos_*'))[0])['consume_benthos']
    #Predator mortality from fishing
    fish_mort_pred = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'fish-mort-pred_*'))[0])['fish_mort_pred']
    #Detritivores mortality from fishing
    fish_mort_det = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'fish-mort-det_*'))[0])['fish_mort_det']
    #Predator mortality from others
    other_mort_pred = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'other-mort-pred_*'))[0])['other_mort_pred']
    #Detritivore mortality from others
    other_mort_det = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'other-mort-det_*'))[0])['other_mort_det']
    #Predator mortality from senescence
    senes_mort_pred = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'senes-mort-pred_*'))[0])['senes_mort_pred']
    #Detritivore mortality from senescence
    senes_mort_det = xr.open_zarr(glob(
        os.path.join(base_folder,
                     'senes-mort-det_*'))[0])['senes_mort_det']

    dbpm_inputs = xr.Dataset(data_vars = {'pref_benthos': pref_benthos,
                                          'pref_pelagic': pref_pelagic,
                                          'constant_growth': constant_growth,
                                          'constant_mortality': constant_mortality,
                                          'consume_pelagic': consume_pelagic,
                                          'consume_benthos': consume_benthos,
                                          'fish_mort_pred': fish_mort_pred,
                                          'fish_mort_det': fish_mort_det,
                                          'other_mort_pred': other_mort_pred,
                                          'other_mort_det': other_mort_det,
                                          'senes_mort_pred': senes_mort_pred,
                                          'senes_mort_det': senes_mort_det})

    return dbpm_inputs


# Loading initialising values for gridded DBPM biomass ------
def loading_dbpm_biomass_inputs(data_folder, init_time = None):
    '''
    Inputs:
    - data_folder (character) Full path to folder where initialising DBPM gridded inputs
    are stored
    - init_time (character) Default is None. Year and month from when to restart gridded 
    DBPM runs. If set to None, DBPM is run from the beginning

    Outputs:
    - ds_init (xarray Dataset) Contains initialising gridded inputs needed to run DBPM
    '''
    if init_time is None:
        predators = xr.open_zarr(
            glob(os.path.join(data_folder, 'predators_*'))[0])['predators'] 
        detritivores = xr.open_zarr(
            glob(os.path.join(data_folder, 'detritivores_*'))[0])['detritivores']
        detritus = xr.open_zarr(
            glob(os.path.join(data_folder, 'detritus_*'))[0])['detritus']
    #If 'init_time' is provided, model runs from 'init_time'
    else:
        predators = xr.open_dataarray(
            glob(os.path.join(data_folder, f'predators_*_{init_time}.nc'))[0])
        detritivores = xr.open_dataarray(
            glob(os.path.join(data_folder, f'detritivores_*_{init_time}.nc'))[0])
        detritus = xr.open_dataarray(
            glob(os.path.join(data_folder, f'detritus_*_{init_time}.nc'))[0])
    
    #Create dataset for predator, detritivores and detritus initialisation data
    ds_init = xr.Dataset(data_vars = {'predators': predators, 
                                      'detritivores': detritivores, 
                                      'detritus': detritus})

    return ds_init


# Loading initialising values for gridded DBPM biomass ------
def loading_dbpm_dynamic_inputs(gridded_esm, gridded_calc, init_time = None, 
                                capped = False):
    '''
    Inputs:
    - gridded_esm (character) Full path to folder where ocean model outputs are stored
    - gridded_calc (character) Full path to folder where processed DBPM inputs are stored
    - init_time (character) Default is None. Year and month from when to restart gridded 
    DBPM runs. If set to None, DBPM is run from the beginning

    Outputs:
    - ds_dynamic (xarray Dataset) Contains dynamic gridded inputs needed to run DBPM
    '''

    if init_time is not None:
        #Get year from initialising time 
        init_yr = pd.Timestamp(init_time).year
        #Timestep from when to restart DBPM 
        subset_time = (pd.Timestamp(init_time)+
                       pd.DateOffset(months = 1)).strftime('%Y-%m')

    if capped:
        cap_search = '-capped'
    else:
        cap_search = ''
    
    if init_time is None or init_yr >= 1841 and init_yr < 1959:
        ui0 = xr.open_mfdataset(glob(os.path.join(
            gridded_calc, f'ui0{cap_search}_spinup*')), engine = 'zarr')['ui0']
        slope = xr.open_mfdataset(glob(os.path.join(
            gridded_esm, f'*spinup_slope{cap_search}_*')), engine = 'zarr')['slope']
        pel_tempeffect = xr.open_mfdataset(glob(
            os.path.join(gridded_calc, 'pel-temp-eff_spinup*')), 
                                           engine = 'zarr')['pel_temp_eff']
        ben_tempeffect = xr.open_mfdataset(glob(os.path.join(
            gridded_calc, 'ben-temp-eff_spinup*')), engine = 'zarr')['ben_temp_eff']
        sinking_rate = xr.open_mfdataset(
            glob(os.path.join(gridded_esm, f'*_spinup_er{cap_search}_*')),
            engine = 'zarr')['export_ratio']
        si_mask = xr.open_mfdataset(glob(os.path.join(
            gridded_esm, f'*_spinup_simask_*')), engine = 'zarr')['simask']
    #Spinup data plus obsclim are loaded if init_time is 1960
    elif init_yr >= 1959 and init_yr < 1961:
        exp = ['spinup', 'obsclim']
        ui0 = xr.open_mfdataset(glob(os.path.join(
            gridded_calc, f'ui0{cap_search}_*')), engine = 'zarr')['ui0']
        slope = xr.open_mfdataset([f for ex in exp for f in glob(os.path.join(
            gridded_esm, f'*{ex}_slope{cap_search}_*'))], engine = 'zarr')['slope']
        pel_tempeffect = xr.open_mfdataset(glob(os.path.join(
            gridded_calc, 'pel-temp-eff_*')), engine = 'zarr')['pel_temp_eff']
        ben_tempeffect = xr.open_mfdataset(glob(os.path.join(
            gridded_calc, 'ben-temp-eff_*')), engine = 'zarr')['ben_temp_eff']
        sinking_rate = xr.open_mfdataset([f for ex in exp for f in glob(os.path.join(
            gridded_esm, f'*{ex}_er{cap_search}_*'))], engine = 'zarr')['export_ratio']
        si_mask = xr.open_mfdataset([f for ex in exp for f in glob(os.path.join(
            gridded_esm, f'*{ex}_simask_*'))], engine = 'zarr')['simask']
    else:
        ui0 = xr.open_mfdataset(glob(os.path.join(
            gridded_calc, f'ui0{cap_search}_[0-9]*')), engine = 'zarr')['ui0']
        slope = xr.open_mfdataset(glob(os.path.join(
            gridded_esm, f'*obsclim_slope{cap_search}_*')), engine = 'zarr')['slope']
        pel_tempeffect = xr.open_mfdataset(glob(os.path.join(
            gridded_calc, 'pel-temp-eff_[0-9]*')), engine = 'zarr')['pel_temp_eff']
        ben_tempeffect = xr.open_mfdataset(glob(os.path.join(
            gridded_calc, 'ben-temp-eff_[0-9]*')), engine = 'zarr')['ben_temp_eff']
        sinking_rate = xr.open_mfdataset(glob(os.path.join(
            gridded_esm, f'*_obsclim_er{cap_search}_*')),
                                         engine = 'zarr')['export_ratio']
        si_mask = xr.open_mfdataset(glob(os.path.join(
            gridded_esm, f'*_obsclim_simask_*')), engine = 'zarr')['simask']

    #Subset data
    if init_time is not None:
        #Subset data from timestep above until the end of the available data
        ui0 = ui0.sel(time = slice(subset_time, None))
        slope = slope.sel(time = slice(subset_time, None))
        pel_tempeffect = pel_tempeffect.sel(time = slice(subset_time, None))
        ben_tempeffect = ben_tempeffect.sel(time = slice(subset_time, None))
        sinking_rate = sinking_rate.sel(time = slice(subset_time, None))
        si_mask = si_mask.sel(time = slice(subset_time, None))

    ds_dynamic = xr.Dataset(data_vars = {'ui0': ui0, 'slope': slope,
                                         'pel_tempeffect': pel_tempeffect,
                                         'ben_tempeffect': ben_tempeffect, 
                                         'sinking_rate': sinking_rate,
                                         'sea_ice_mask': si_mask})

    return ds_dynamic


# Calculating feeding and satiation rates ------
def feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                            dbpm_dynamic_inputs, group):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    - dbpm_init_inputs (xarray Dataset) Contains gridded inputs with initialisation 
    values for predators, detritivores and detritus
    - dbpm_dynamic_inputs (xarray Dataset) Contains dynamic gridded inputs need to run
    DBPM
    - group (character): Choose the group for which feeding and satiation rates will 
    be calculated for. Choices: predators, detritivores, detritus.

    Outputs:
    - feed_sat_rates (xarray Dataset). Contains feeding and satiation rates for chosen
    group
    '''
    
    if group == 'predators':
        #To be applied to feeding rates for pelagics and benthic groups
        feed_mult_pel = (gridded_params['hr_volume_search']*
                         (10**(dbpm_fixed_inputs['log10_size_bins']*
                               gridded_params['metabolic_req_pred']))*
                         dbpm_fixed_inputs['pref_pelagic'])
        # Predators (f_pel)
        pred_growth = (((dbpm_init_inputs['predators']*
                         gridded_params['log_size_increase']).
                        dot(dbpm_fixed_inputs['constant_growth']).
                       rename({'sc': 'size_class'}))*feed_mult_pel)
        feed_rate_pel = (dbpm_dynamic_inputs['pel_tempeffect']*
                         (pred_growth/(1+gridded_params['handling']*pred_growth)))
        # Satiation level of predator for pelagic prey (sat.pel)
        sat_pel = xr.where(feed_rate_pel > 0, feed_rate_pel/pred_growth, 0)

        #Create dataset with feeding and satiation outputs for predators
        feed_sat_rates = xr.Dataset(data_vars = {'f_pel': feed_rate_pel,
                                                 'sat_pel': sat_pel})
    
    elif group == 'detritivores':
        #To be applied to feeding rates for benthic group
        feed_mult_ben = (gridded_params['hr_volume_search']*
                         (10**(dbpm_fixed_inputs['log10_size_bins']*
                               gridded_params['metabolic_req_pred']))*
                        dbpm_fixed_inputs['pref_benthos'])
        # Detritivores (f_ben)
        detrit_growth = (((dbpm_init_inputs['detritivores']*
                           gridded_params['log_size_increase']).
                          dot(dbpm_fixed_inputs['constant_growth']).
                          rename({'sc': 'size_class'}))*feed_mult_ben)
        feed_rate_bent = (dbpm_dynamic_inputs['pel_tempeffect']*
                          (detrit_growth/(1+gridded_params['handling']*detrit_growth)))
        # Satiation level of predator for benthic prey
        calc_growth = ((gridded_params['hr_volume_search']*
                        (10**(dbpm_fixed_inputs['log10_size_bins']*
                              gridded_params['metabolic_req_detritivore']))*
                        dbpm_fixed_inputs['pref_benthos'])*
                       ((dbpm_init_inputs['detritivores']*
                         gridded_params['log_size_increase']).
                        dot(dbpm_fixed_inputs['constant_growth'])).
                       rename({'sc': 'size_class'}))
        #Satiation level of predator for benthic prey (sat.ben)
        sat_ben = xr.where(feed_rate_bent > 0, feed_rate_bent/calc_growth, 0)

        #Create dataset with feeding and satiation outputs for detritivores
        feed_sat_rates = xr.Dataset(data_vars = {'f_ben': feed_rate_bent,
                                                 'sat_ben': sat_ben})

    elif group == 'detritus':
        # Feeding rates
        detritus_multiplier = ((1/dbpm_fixed_inputs['size_bin_vals'])*
                               gridded_params['hr_vol_filter_benthos']*
                               10**(dbpm_fixed_inputs['log10_size_bins']*
                                    gridded_params['metabolic_req_detritivore'])*
                               dbpm_init_inputs['detritus'])
        feed_rate_det = (dbpm_dynamic_inputs['ben_tempeffect']*detritus_multiplier/
                         (1+gridded_params['handling']*detritus_multiplier))
        
        #Create dataset with feeding and satiation outputs
        feed_sat_rates = xr.Dataset(data_vars = {'f_det': feed_rate_det})

    else:
        print('the "group" parameter must be "predators", "detritivores", or "detritus"')

    #Reorganise output dimensions
    feed_sat_rates = feed_sat_rates.transpose('time', 'size_class', 'lat', 'lon')
    #Apply spatial mask
    feed_sat_rates = feed_sat_rates.where(dbpm_fixed_inputs['mask'])
    #Ensure correct time is applied
    feed_sat_rates['time'] = dbpm_init_inputs.time
    
    return feed_sat_rates


def mortality_calc(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                   dbpm_dynamic_inputs, feed_sat_rates, group):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    - dbpm_init_inputs (xarray Dataset) Contains gridded inputs with initialisation 
    values for predators, detritivores and detritus
    - dbpm_dynamic_inputs (xarray Dataset) Contains dynamic gridded inputs need to run
    DBPM
    - feed_sat_rates (xarray Dataset) Contains feeding and satiation rates for the 
    group for which mortality rates will be calculated. This is the output of the
    "feeding_satiation_rates" function
    - group (character): Choose the group for which mortality rates will be calculated.
    Choices: predators, detritivores.

    Outputs:
    - mortality_terms (xarray Dataset). Contains mortality terms for all groups modelled
    '''

    if group == 'predators': 
        # Fishing mortality (FVec.u, FVec.v) from Benoit & Rochet 2004. Here fish_mort_pred 
        # and fish_mort_pred= fixed catchability term for predators and detritivores to be 
        # estimated along with ind_min_det and ind_min_fish_pred
        fishing_mort_pred = (dbpm_fixed_inputs['fish_mort_pred']*
                             dbpm_dynamic_inputs['effort']).drop_vars('time')
        #Time dimension matches feed rates time stamp
        fishing_mort_pred = fishing_mort_pred.expand_dims({'time': dbpm_init_inputs.time})
        # Predator death integrals
        # Predation mortality
        # Predators (PM.u)
        pred_mort_pred = (dbpm_fixed_inputs['consume_pelagic']*
                          ((dbpm_init_inputs['predators']*feed_sat_rates['sat_pel']*
                             gridded_params['log_size_increase']).
                           dot(dbpm_fixed_inputs['constant_mortality'])).
                          rename({'sc': 'size_class'}))
        # Total mortality
        # Predators (Z.u)
        tot_mort_pred = (pred_mort_pred+dbpm_dynamic_inputs['pel_tempeffect']*
                         dbpm_fixed_inputs['other_mort_pred']+
                         dbpm_fixed_inputs['senes_mort_pred']+fishing_mort_pred)
        
        # Create dataset with mortality for predators
        mortality_terms = xr.Dataset(data_vars = {'Fvec_u': fishing_mort_pred,
                                                  'PM_u': pred_mort_pred,
                                                  'Z_u': tot_mort_pred})
    elif group == 'detritivores':
        # Fishing mortality (FVec.u, FVec.v) from Benoit & Rochet 2004. Here 
        # fish_mort_pred and fish_mort_pred= fixed catchability term for predators and 
        # detritivores to be estimated along with ind_min_det and ind_min_fish_pred
        fishing_mort_det = (dbpm_fixed_inputs['fish_mort_det']*
                            dbpm_dynamic_inputs['effort']).drop_vars('time')
        #Time dimension matches feed rates time stamp
        fishing_mort_det = fishing_mort_det.expand_dims({'time': dbpm_init_inputs.time})
        # Predation mortality
        # Detritivores (PM.v)
        pred_mort_det = xr.where(feed_sat_rates['sat_ben'] > 0, 
                                 (dbpm_fixed_inputs['consume_benthos']*
                                  ((dbpm_init_inputs['predators']*
                                    feed_sat_rates['sat_ben']*
                                    gridded_params['log_size_increase']).
                                   dot(dbpm_fixed_inputs['constant_mortality'])).
                                  rename({'sc': 'size_class'})), 0)
        # Total mortality
        # Detritivores (Z.v)
        tot_mort_det = (pred_mort_det+dbpm_dynamic_inputs['ben_tempeffect']*
                        dbpm_fixed_inputs['other_mort_det']+
                        dbpm_fixed_inputs['senes_mort_det']+fishing_mort_det)

        # Create dataset with mortality for detritivores
        mortality_terms = xr.Dataset(data_vars = {'Fvec_v': fishing_mort_det,
                                                  'PM_v': pred_mort_det,
                                                  'Z_v': tot_mort_det})
    else:
        print('the "group" parameter must be either "predators" or "detritivores"')
        
    #Reorganise dimensions
    mortality_terms = mortality_terms.transpose('time', 'size_class', 'lat', 'lon')
    #Apply spatial mask
    mortality_terms = mortality_terms.where(dbpm_fixed_inputs['mask'])
    #Ensure correct time is applied
    mortality_terms['time'] = dbpm_init_inputs.time
    
    return mortality_terms


def growth_rates_calc(gridded_params, feed_sat_rates, dbpm_fixed_inputs):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - feed_sat_rates (xarray Dataset) Contains feeding rates calculated for predators,
    detritivores and detritus calculated with the "feeding_satiation_rates" function
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    
    Outputs:
    - growth_rates (xarray Dataset). Contains growth and estimates for detritivores 
    and predators
    '''
    # Predator growth integral (GG_u)
    growth_int_pred = (gridded_params['growth_prop']*gridded_params['growth_pred']*
                       feed_sat_rates['f_pel']+gridded_params['high_prop']*
                       gridded_params['growth_detritivore']*feed_sat_rates['f_ben'])
    # Benthos growth integral yr-1 (GG_v)
    growth_int_det = (gridded_params['high_prop']*gridded_params['growth_detritus']*
                      feed_sat_rates['f_det'])
    # Create dataset with growth rates for detritivores and predators
    growth_reprod = xr.Dataset(data_vars = {'GG_u': growth_int_pred,
                                            'GG_v': growth_int_det})

    #Reorganise dimensions
    growth_reprod = growth_reprod.transpose('time', 'size_class', 'lat', 'lon')
    #Apply spatial mask
    growth_reprod = growth_reprod.where(dbpm_fixed_inputs['mask'])
    #Ensure correct time is applied
    growth_reprod['time'] = feed_sat_rates.time
    
    return growth_reprod


def reproduction_rates(gridded_params, feed_sat_rates, dbpm_fixed_inputs, group):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - feed_sat_rates (xarray Dataset) Contains feeding rates calculated for predators,
    detritivores and detritus calculated with the "feeding_satiation_rates" function
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    - group (character): Choose the group for which mortality rates will be 
    calculated. Choices: predators, detritivores.
    
    Outputs:
    - growth_rates (xarray Dataset). Contains growth and estimates for detritivores 
    and predators
    '''

    if group == 'predators':
        #Keep relevant size classes
        feed_sat_rates = (feed_sat_rates.
                          isel(size_class = 
                               slice(gridded_params['ind_min_pred_size']+1, None)))
        #Predators
        reprod_pred = (gridded_params['growth_prop']*
                       (1-(gridded_params['growth_pred']+
                           gridded_params['energy_pred']))*
                       feed_sat_rates['f_pel']+gridded_params['growth_prop']*
                       (1-(gridded_params['growth_detritivore']+
                           gridded_params['energy_detritivore']))*
                       feed_sat_rates['f_ben'])
        #Adding to output dataset
        reprod_rates = xr.Dataset(data_vars = {'R_u': reprod_pred})
        
    elif group == 'detritivores':
        #Keep relevant size classes
        feed_sat_rates = (feed_sat_rates.
                          isel(size_class = 
                               slice(gridded_params['ind_min_detritivore_size']+1, 
                                     None)))
        #Detritivores
        reprod_det = (gridded_params['high_prop']*
                      (1-(gridded_params['growth_detritus']+
                          gridded_params['energy_detritivore']))*
                      feed_sat_rates['f_det'])
        #Adding to output dataset
        reprod_rates = xr.Dataset(data_vars = {'R_v': reprod_det})

    else:
        print('the "group" parameter must be either "predators" or "detritivores"')

    #Reorganise dimensions
    reprod_rates = reprod_rates.transpose('time', 'size_class', 'lat', 'lon')
    #Apply spatial mask
    reprod_rates = reprod_rates.where(dbpm_fixed_inputs['mask'])
    #Ensure correct time is applied
    reprod_rates['time'] = feed_sat_rates.time
    
    return reprod_rates


def detritus_output(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                    feed_detritus):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    - dbpm_init_inputs (xarray Dataset) Contains gridded inputs with initialisation 
    values for predators, detritivores and detritus
    - feed_detritus (xarray Data Array) Contains detritus rates calculated with the 
    "feeding_satiation_rates" function
    
    Outputs:
    - output_w (xarray Data Array). Contains detritus output (g.m-2.yr-1), which are
    defined as losses from detritivore scavenging/filtering only.
    '''
    # Detritus output (g.m-2.yr-1)
    # losses from detritivore scavenging/filtering only:
    output_w = (dbpm_fixed_inputs['size_bin_vals']*
                gridded_params['log_size_increase']*
                dbpm_init_inputs['detritivores']*
                feed_detritus).sum('size_class')

    # Adding variable name
    output_w.name = 'output_w'
    #Reorganise dimensions
    output_w = output_w.transpose('time', 'lat', 'lon')
    #Apply spatial mask
    output_w = output_w.where(dbpm_fixed_inputs['mask'])
    
    return output_w
    

def defecation_predators(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                         feed_sat_rates):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    - dbpm_init_inputs (xarray Dataset) Contains gridded inputs with initialisation 
    values for predators, detritivores and detritus
    - feed_sat_rates (xarray Dataset) Contains feeding rates for predators and 
    detritivores calculated with the "feeding_satiation_rates" function
    
    Outputs:
    - defbypred (xarray Data Array). Contains total biomass defecated by predators
    '''
    # Select relevant predator size classes
    pred = (dbpm_init_inputs['predators'].
            isel(size_class = slice(gridded_params['ind_min_pred_size'], None)))

    # Total biomass density defecated by pred (g.m-2.yr-1)
    defbypred = ((gridded_params['defecate_prop']*feed_sat_rates['f_pel']*
                  dbpm_fixed_inputs['size_bin_vals']*
                  pred+gridded_params['def_low']*feed_sat_rates['f_ben']*
                  dbpm_fixed_inputs['size_bin_vals']*pred)*
                 gridded_params['log_size_increase']).sum('size_class')  

    # Adding variable name
    defbypred.name = 'defbypred'
    #Reorganise dimensions
    defbypred = defbypred.transpose('time', 'lat', 'lon')
    #Apply spatial mask
    defbypred = defbypred.where(dbpm_fixed_inputs['mask'])
    
    return defbypred


def detritus_pool(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                  dbpm_dynamic_inputs, defbypred, output_w):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to 
    run DBPM
    - dbpm_init_inputs (xarray Dataset) Contains gridded inputs with initialisation 
    values for predators, detritivores and detritus
    - dbpm_dynamic_inputs (xarray Dataset) Contains dynamic gridded inputs need to
    run DBPM
     - defbypred (xarray Data Array). Contains total biomass defecated by predators 
     calculated with the "defecation_predators" function
     - output_w (xarray Data Array). Contains detritus output (g.m-2.yr-1) 
     calculated with the "detritus_output" function
    
    Outputs:
    - det_pool (xarray Dataset). Contains detritus biomass density pool fluxes
    (g.m-2.yr-1)
    '''

    #Mltipliers used more than once in the input_w calculation
    det_multiplier = (dbpm_dynamic_inputs['ben_tempeffect']*
                      dbpm_fixed_inputs['size_bin_vals']*
                      gridded_params['log_size_increase']*
                      dbpm_init_inputs['detritivores'])

    # Considering pelagic faeces as input as well as dead bodies 
    # from both pelagic and benthic communities and phytodetritus
    # (dying sinking phytoplankton)
    if gridded_params['detritus_coupling']:
        pred_multiplier = (dbpm_dynamic_inputs['pel_tempeffect']*
                           dbpm_init_inputs['predators']*
                           dbpm_fixed_inputs['size_bin_vals']*
                           gridded_params['log_size_increase'])
        # Pelagic spectrum inputs (sinking dead bodies and faeces): 
        # Export ratio used for "sinking rate" + benthic spectrum 
        # inputs (dead stuff already on/in seafloor)
        input_w = (dbpm_dynamic_inputs['sinking_rate']*
                   (defbypred+
                    ((dbpm_fixed_inputs['other_mort_pred']*pred_multiplier)+ 
                     (dbpm_fixed_inputs['senes_mort_pred']*pred_multiplier)+
                     (dbpm_fixed_inputs['other_mort_det']*det_multiplier)+
                     (dbpm_fixed_inputs['senes_mort_det']*det_multiplier))
                    .sum('size_class')))
    else:
        input_w = (((dbpm_fixed_inputs['other_mort_det']*det_multiplier)+ 
                    (dbpm_fixed_inputs['senes_mort_det']*det_multiplier)).
                   sum('size_class'))
   
   # Get burial rate from Dunne et al. 2007 equation 3
    burial = input_w*(0.013+0.53*input_w**2/(7+input_w)**2)
    
    # Losses from detritivory + burial rate (not including 
    # remineralisation because that goes to p.p. after sediment, 
    # we are using realised p.p. as inputs to the model) 
    dW = input_w-(output_w+burial)
   
    # Create dataset with changes in detritus biomass
    det_pool = xr.Dataset(data_vars = {'input_w': input_w, 
                                       'burial': burial, 
                                       'dW': dW})
    #Reorganise dimensions
    det_pool = det_pool.transpose('time', 'lat', 'lon')
    #Apply spatial mask
    det_pool = det_pool.where(dbpm_fixed_inputs['mask'])
    
    return det_pool


def biomass_density(gridded_params, dbpm_fixed_inputs, growth_rate, 
                     total_mortality, biomass, group):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to 
    run DBPM
    - growth_rates (xarray Dataset). Contains growth and estimates for detritivores 
    and predators estimated by the "growth_rate" function
    - total_mortality (xarray Data Array). Contains total mortality for group of
    interest estimated by the "mortality_calc" function
    - biomass (xarray Data Array). Contains total biomass for group of interest
    - group (character): Choose the group for which mortality rates will be 
    calculated. Choices: predators, detritivores.
    
    Outputs:
    - density (xarray Dataset). Contains biomass density (nos.m-2) for the group
    of interest
    '''
    
    if group == 'predators':
        index_sc = np.array(gridded_params['idx'])
        ref_sc =  gridded_params['log10_ind_min_pred_size']
    elif group == 'detritivores':
        index_sc = np.array(gridded_params['idx_new'])
        ref_sc =  gridded_params['log10_ind_min_detritivore_size']
    else:
        print('the "group" parameter must be either "predators" or ' + 
              '"detritivores"')
    
    # Pelagic Predator Density (nos.m-2) ----
    # Solve for time + delta.t using implicit time Euler upwind finite 
    # difference (help from Ken Andersen and Richard Law)

    # Shifting size classes by one before calculating Ai
    ggp_shift = growth_rate.shift({'size_class': 1})
    ggp_shift, range_sc = xr.align(ggp_shift, 
                                   dbpm_fixed_inputs['log10_size_bins'].
                                   isel(size_class = index_sc))
      
    # Matrix setup for implicit differencing 
    Ai = (((1/np.log(10))*-ggp_shift*gridded_params['timesteps_years']/
           gridded_params['log_size_increase']).
          reindex(size_class = dbpm_fixed_inputs['log10_size_bins'], 
                  fill_value = 0))     
    
    Bi = ((1+(1/np.log(10))*growth_rate.isel(size_class = index_sc)*
           gridded_params['timesteps_years']/
           gridded_params['log_size_increase']+
           total_mortality.isel(size_class = index_sc)*
             gridded_params['timesteps_years']).
            reindex(size_class = dbpm_fixed_inputs['log10_size_bins'],
                    fill_value = 0))
    
    Si = (biomass.isel(size_class = index_sc).
          reindex(size_class = dbpm_fixed_inputs['log10_size_bins'],
                  fill_value = 0))

    # Boundary condition at upstream end
    Ai = xr.where(Ai.size_class == ref_sc, 0, Ai)   
    Bi = xr.where(Bi.size_class == ref_sc, 1, Bi)
    if ref_sc not in growth_rate.size_class[index_sc]:
        Si = xr.where(Si.size_class == ref_sc, biomass, Si)

    density = xr.Dataset(data_vars = {'Ai': Ai, 'Bi': Bi, 'Si': Si})
    #Reorganise dimensions
    density = density.transpose('time', 'size_class', 'lat', 'lon')
    #Apply spatial mask
    density = density.where(dbpm_fixed_inputs['mask'])

    return density


def tot_biomass_calc(gridded_params, dbpm_fixed_inputs, group, biomass_current,
                     biomass_next, biomass_density, reprod_rate = None, 
                     growth_rate = None, total_mortality = None):

    if group == 'predators':
        ref_sc =  gridded_params['log10_ind_min_pred_size']
        ref = gridded_params['ind_min_pred_size']
    elif group == 'detritivores':
        ref_sc =  gridded_params['log10_ind_min_detritivore_size']
        ref = gridded_params['ind_min_detritivore_size']
    else:
        print('the "group" parameter must be either "predators" or ' +
              '"detritivores"')

    # Apply transfer efficiency of 10% *plankton density at same size
    # reproduction from energy allocation
    if gridded_params['dynamic_reproduction']:
        if growth_rate is None:
            raise Exception('When "dynamic_reproduction" is set to True, ' + 
                            '"growth_rate" must be provided as input')
        elif total_mortality is None:
            raise Exception('When "dynamic_reproduction" is set to True, ' +
                            '"total_mortality" must be provided as input')
        elif reprod_rate is None:
            raise Exception('When "dynamic_reproduction" is set to True, ' + 
                            '"reprod_rate" must be provided as input')
        else:
            growth_rate = growth_rate.sel(size_class = ref_sc)
            total_mortality = total_mortality.sel(size_class = ref_sc)
        bio_ref = biomass_current.sel(size_class = ref_sc)
        bio_ref_next = (bio_ref+
                        ((reprod_rate*dbpm_fixed_inputs['size_bin_vals']*
                          biomass_current*gridded_params['log_size_increase']).
                         sum('size_class')*gridded_params['timesteps_years'])/
                        (gridded_params['log_size_increase']*
                         dbpm_fixed_inputs['size_bin_vals'].sel(size_class = ref_sc))-
                        (gridded_params['timesteps_years']/
                         gridded_params['log_size_increase'])*(1/np.log(10))*
                        growth_rate*bio_ref-gridded_params['timesteps_years']*
                        total_mortality*bio_ref)
        # Change time stamp to next month
        if isinstance(biomass_next.time.values, np.ndarray):
            bio_ref_next['time'] = biomass_next.time.values
        else:
            bio_ref_next['time'] = [biomass_next.time.values]
        #Replace boundary condition in next timestep
        biomass_next = xr.where(biomass_next.size_class == ref_sc, bio_ref_next,
                                biomass_next)
        del bio_ref, bio_ref_next

    #main loop calculation
    biomass_density = biomass_density.drop_vars('time').squeeze()
    for sc in range(ref+1, gridded_params['numb_size_bins']):
        bio_next_shift = biomass_next.isel(size_class = sc-1).drop_vars('size_class')
        da = ((biomass_density['Si'].isel(size_class = sc)-
               biomass_density['Ai'].isel(size_class = sc)*bio_next_shift)/
              biomass_density['Bi'].isel(size_class = sc))
        biomass_next = xr.where(biomass_next.size_class == da.size_class, da, 
                                biomass_next)
    
    #Adding name to data array
    biomass_next.name = group
    #Reorganise dimensions
    biomass_next = biomass_next.transpose('time', 'size_class', 'lat', 'lon')
    #Apply spatial mask
    biomass_next = biomass_next.where(dbpm_fixed_inputs['mask'])
        
    return biomass_next
    

# Run model per grid cell or averaged over an area ------
def gridded_sizemodel(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs, 
                      dbpm_dynamic_inputs, region, model_res, out_folder):
    '''
    Inputs:
    - gridded_params (dictionary) Gridded DBPM parameters obtained in step 04.
    - dbpm_fixed_inputs (xarray Dataset) Contains fixed gridded inputs needed to run 
    DBPM
    - dbpm_init_inputs (xarray Dataset) Contains gridded inputs with initialisation 
    values for predators, detritivores and detritus
    - dbpm_dynamic_inputs (xarray Dataset) Contains dynamic gridded inputs need to run
    DBPM
    - region (character) Name of region being processed (as included in folder and file
    names)
    - model_res (character) Spatial resolution of DBPM inputs (as included in folder
    and file names)
    - out_folder (character) Full path to folder where DBPM outputs will be stored

    Outputs:
    - None. This function does not return any objects. Instead outputs are saved in
    the path provided in the "out_folder" argument
    '''
    
    # Storing date from dynamic dataset
    dbpm_time = dbpm_dynamic_inputs.time.values
    pred_ts_next = pd.Timestamp(dbpm_time).strftime('%Y-%m')

    # Storing date from initialising dataset
    [pred_time] = dbpm_init_inputs.time.values
    pred_ts = pd.Timestamp(pred_time).strftime('%Y-%m')

    # Feeding and satiation rates ----
    # Predators
    feed_sat_pred = feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, 
                                            dbpm_init_inputs, dbpm_dynamic_inputs, 
                                            group = 'predators')
    # Detritivores
    feed_sat_det = feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, 
                                           dbpm_init_inputs, dbpm_dynamic_inputs, 
                                           group = 'detritivores')
    # Detritus
    feed_detritus = feeding_satiation_rates(gridded_params, dbpm_fixed_inputs, 
                                            dbpm_init_inputs, dbpm_dynamic_inputs, 
                                            group = 'detritus')
    
    
    # Mortality terms ----
    # Predators
    mortality_pred = mortality_calc(gridded_params, dbpm_fixed_inputs, 
                                    dbpm_init_inputs, dbpm_dynamic_inputs, 
                                    feed_sat_pred, group = 'predators')
    # Saving mortality rates for predators
    fn = f'pred_mort_pred_{model_res}_{region}_{pred_ts}.nc'
    mortality_pred['PM_u'].to_netcdf(os.path.join(out_folder, fn))

    # Detritivores
    mortality_det = mortality_calc(gridded_params, dbpm_fixed_inputs, 
                                   dbpm_init_inputs, dbpm_dynamic_inputs, 
                                   feed_sat_det, group = 'detritivores')
    # Saving mortality rates for detritivores
    fn = f'pred_mort_det_{model_res}_{region}_{pred_ts}.nc'
    mortality_det['PM_v'].to_netcdf(os.path.join(out_folder, fn))

    
    # Growth rates ----
    # Predators and detritivores
    growth_rates_pred_det = growth_rates_calc(gridded_params, 
                                              xr.Dataset({'f_pel': feed_sat_pred.f_pel,
                                                          'f_ben': feed_sat_det.f_ben,
                                                          'f_det': feed_detritus.f_det}), 
                                              dbpm_fixed_inputs)
    # Saving growth rates for predators
    fn = f'growth_int_pred_{model_res}_{region}_{pred_ts}.nc'
    growth_rates_pred_det['GG_u'].to_netcdf(os.path.join(out_folder, fn))
    # Saving growth rates for detritivores
    fn = f'growth_int_det_{model_res}_{region}_{pred_ts}.nc'
    growth_rates_pred_det['GG_v'].to_netcdf(os.path.join(out_folder, fn))


    # Reproduction rates ----
    if gridded_params['dynamic_reproduction']:
        reprod_pred = reproduction_rates(gridded_params, 
                                         xr.Dataset({'f_pel': feed_sat_pred.f_pel,
                                                     'f_ben': feed_sat_det.f_ben}),
                                         dbpm_fixed_inputs, group = 'predators')
        reprod_det = reproduction_rates(gridded_params, feed_detritus,
                                        dbpm_fixed_inputs, group = 'detritivores')
    # Detritus output (g.m-2.yr-1) ----
    # losses from detritivore scavenging/filtering only:
    output_w = detritus_output(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs,
                               feed_detritus.f_det)


    # Defecation by predators ----
    defbypred = defecation_predators(gridded_params, dbpm_fixed_inputs,
                                     dbpm_init_inputs, 
                                     xr.Dataset({'f_pel': feed_sat_pred.f_pel,
                                                 'f_ben': feed_sat_det.f_ben}))

    # Detritus Biomass Density Pool ----
    # Fluxes in and out (g.m-2.yr-1) of detritus pool. Solve for detritus biomass
    # density in next timestep 
    # Increment values of detritus, predators & detritivores for next 
    # timestep
    det_pool = detritus_pool(gridded_params, dbpm_fixed_inputs, dbpm_init_inputs,
                             dbpm_dynamic_inputs, defbypred, output_w)

    # Biomass density of detritus g.m-2
    detritus = (dbpm_init_inputs['detritus']+det_pool['dW']*
                gridded_params['timesteps_years']).load()

    # Updating timestamp (results for next time step)
    detritus['time'] = [dbpm_time]
    #Adding name
    detritus.name = 'detritus'
    # Saving detritus biomass
    fn = f'detritus_{model_res}_{region}_{pred_ts_next}.nc'
    detritus.to_netcdf(os.path.join(out_folder, fn))

    
    # Pelagic predator density ----
    pred_density = biomass_density(gridded_params, dbpm_fixed_inputs, 
                                   growth_rates_pred_det['GG_u'],
                                   mortality_pred['Z_u'], 
                                   dbpm_init_inputs['predators'], 
                                   group = 'predators')

    
    # Recruitment at smallest consumer mass ----
    # continuation of plankton hold constant  
    pred_next = (dbpm_dynamic_inputs['ui0']*
                 (10**(dbpm_dynamic_inputs['slope']*
                       dbpm_fixed_inputs['log10_size_bins'])))
    pred_next = xr.where((pred_next.size_class <= 
                          gridded_params['log10_ind_min_pred_size']),
                         pred_next, 0)

    predators = tot_biomass_calc(gridded_params, dbpm_fixed_inputs, 
                                 'predators', dbpm_init_inputs['predators'], pred_next, 
                                 pred_density, reprod_rate = reprod_pred['R_u'],
                                 growth_rate = growth_rates_pred_det['GG_u'], 
                                 total_mortality = mortality_pred['Z_u']).load()
    #If values are negative, assign a value of 0
    predators = xr.where(predators < 0, 0, predators)
    
    # Saving predators biomass
    # Creating file name
    fn = f'predators_{model_res}_{region}_{pred_ts_next}.nc'
    predators.to_netcdf(os.path.join(out_folder, fn))

    
    # Benthic detritivore density ----
    det_density = biomass_density(gridded_params, dbpm_fixed_inputs, 
                                  growth_rates_pred_det['GG_v'], 
                                  mortality_det['Z_v'], 
                                  dbpm_init_inputs['detritivores'], 
                                  group = 'detritivores')

    
    # Recruitment at smallest detritivore mass ----
    # hold constant continution of plankton with sinking rate multiplier 
    detriti_next = xr.where((dbpm_init_inputs['detritivores'].size_class <= 
                             gridded_params['log10_ind_min_detritivore_size']),
                            dbpm_init_inputs['detritivores'], 0)
    # Updating timestamp (results for next time step)
    detriti_next['time'] = [dbpm_time]

    detritivores = tot_biomass_calc(gridded_params, dbpm_fixed_inputs, 
                                    'detritivores', dbpm_init_inputs['detritivores'],
                                    detriti_next, det_density, 
                                    reprod_rate = reprod_det['R_v'],
                                    growth_rate = growth_rates_pred_det['GG_v'], 
                                    total_mortality = mortality_det['Z_v']).load()
    #If values are negative, assign a value of 0
    detritivores = xr.where(detritivores < 0, 0, detritivores)
    
    # Saving detritivores biomass
    fn = f'detritivores_{model_res}_{region}_{pred_ts_next}.nc'
    detritivores.to_netcdf(os.path.join(out_folder, fn))

    #Creating dataset with new values for predators and detritivores biomass, and detritus
    ds_next = xr.Dataset(data_vars = {'predators': predators, 
                                      'detritivores': detritivores, 
                                      'detritus': detritus})

    return ds_next


# Calculate decadal catches ------
def decadal_catches(catches_da, mask_da, **kwargs):
    '''
    Inputs:
    - catches_da (data array) 3D data array containing total catches over time. 
    - mask_da (data array) 2D data array used to mask out any areas outside the
    region of interest. Areas outside boundaries must be identified with NaN.
    **Optional**: 
    - name (character) Name to be applied to decadal data array
    - attrs (dictionary) Containing attributes to be added to data array

    Output:
    - dec_da (data array): Mean catches per decade
    '''

    catches_da['decade'] = (np.floor(catches_da.time.dt.year/10)*10).astype(int)
    dec_da = catches_da.groupby('decade').mean()
    dec_da = dec_da.where(np.isfinite(mask_da))
    if 'name' in kwargs.keys():
        dec_da.name = kwargs.get('name')
    if 'attrs' in kwargs.keys():
        dec_da = dec_da.assign_attrs(kwargs.get('attrs'))
    return dec_da


# Calculate lowest maximum difference across timesteps ----
def lowest_maximum_absolute_diff(folder_path, resolution, var_name, prop = 1):
    '''
    folder_path (character) - Full path to folder containing files of interest
    resolution (character) - Resolution of the data to be used in calculations
    var_name (character) - Name of variable of interest as included in file name
    prop (numeric) - Default is 1. Multiplier to be applied to lowest maximum
    absolute difference if problems with model persist

    Output:
    low_max_diff (numeric) - Smallest maximum absolute difference across time
    '''
    # Find ctrlclim and obsclim files
    [ctrl_file] = glob(os.path.join(folder_path, resolution, f'*_ctrlclim_{var_name}_*'))
    [obs_file] = glob(os.path.join(folder_path, resolution, f'*_obsclim_{var_name}_*'))
    # Load files
    ds_ctrl = abs(xr.open_zarr(ctrl_file)[var_name].diff(dim = 'time'))
    ds_obs = abs(xr.open_zarr(obs_file)[var_name].diff(dim = 'time'))
    # Calculate lowest maximum difference
    low_max_diff = min(ds_ctrl.max().values, ds_obs.max().values)
    return low_max_diff*prop


# Remove grid cells exceeding temporal difference
def mask_absolute_diff(folder_path, resolution, experiment, var_name, threshold):
    '''
    folder_path (character) - Full path to folder containing files of interest
    resolution (character) - Resolution of the data to be used in calculations
    experiment (character) - Experiment name: obsclim, ctrlclim, spinup, or 
    stable-spin
    var_name (character) - Name of variable of interest as included in file name
    threshold (numeric) - Threshold is used as the maximum temporal difference 
    allowed in the dataset

    Output:
    No outputs returned. Masked xarray dataset saved to the same folder provided  
    in 'folder_path'.
    '''
    # Identify file name
    [file_target] = glob(os.path.join(folder_path, resolution, f'*_{experiment}_{var_name}_*'))
    # Creating filename to store result
    fout = os.path.basename(file_target).replace(f'_{var_name}_', f'_{var_name}-capped_')
    # Load data
    ds = xr.open_zarr(file_target)[var_name]
    # Calculate absolute difference across time steps
    ds_diff = abs(ds.diff(dim = 'time'))
    # Identify grid cells where threshold is exceeded
    mask = xr.where(ds > threshold, 1, 0)
    # Assign NA to grid cells where threshold is exceeded at least once
    ds_masked = xr.where(mask > 0, np.nan, ds)
    # Save result
    fout = os.path.join(folder_path, resolution, fout)
    ds_masked.to_zarr(fout, consolidated = True, mode = 'w')
    
    numb_nas = mask.sum().values
    print(f'Total number of grid cells exceeding threshold ({threshold}) '+
          f'in {var_name} {experiment}: {numb_nas}')


# Recalculating phytoplankton-based variables from masked phytoplankton data 
def exportRatio_intercept_slope(folder_path, resolution, experiment,
                                mmin = 10**(-14.25), mmid = 10**(-10.184), 
                                mmax = 10**(-5.25)):
    '''
    Inputs:
    - folder_path (character) File path pointing to folder containing
    zarr files with GFDL data for the region of interest
    - resolution (character) Resolution of the data that will be processed
    - experiment (character) Select 'ctrlclim','obsclim', 'spinup', 'stable-spin'
    The parameters below come from the 'GetPPIntSlope' function:
    - mmin (numeric)  Default is 10**(-14.25). ????
    - mmid (numeric)  Default is 10**(-10.184). ????
    - mmax (numeric)  Default is 10**(-5.25). ????

    Outputs:
    No outputs returned, but export ratio, as well as the slope and intercept for 
    phytoplankton are saved in the folder provided in 'folder_path'. 
    ui0 is also calculated and saved, but under 'gridded_params'
    '''

    #Get list of files in experiment
    file_list = glob(os.path.join(folder_path, resolution, f'*_{experiment}_*'))

    #load depth
    if experiment == 'obsclim':
        [depth_file] = glob(os.path.join(folder_path, resolution, 
                                         f'*_{experiment}_deptho_*'))
    else:
        [depth_file] = glob(os.path.join(folder_path, resolution,
                                         f'*_ctrlclim_deptho_*'))
    depth = xr.open_zarr(depth_file)['deptho']
    
    #Load sea surface temperature
    tos = xr.open_zarr([f for f in file_list if '_tos_' in f][0])['tos']
    
    #Load small phytoplankton
    sphy = xr.open_zarr([f for f in file_list if '_sphy-capped_' in f][0])['sphy']
    
    #Load large phytoplankton
    lphy = xr.open_zarr([f for f in file_list if '_lphy-capped_' in f][0])['lphy']

    #Calculate total phytoplankton
    ptotal = lphy+sphy

    #Calculate phytoplankton size ratios
    plarge = lphy/ptotal
    psmall = sphy/ptotal

    #Calculate export ration
    er = (np.exp(-0.032*tos)*((0.14*psmall)+(0.74*(plarge)))+
          (0.0228*(plarge)*(depth*0.004)))/(1+(depth*0.004))
    #If values are negative, assign a value of 0
    er = xr.where(er < 0, 0, er)
    #If values are above 1, assign a value of 1
    er = xr.where(er > 1, 1, er)
    er.name = 'export_ratio'

    #Creating file path to save export ratio
    [fout] = [os.path.basename(f) for f in file_list if '_er_' in f]
    fout = os.path.join(folder_path, '025deg', fout.replace('_er_', 
                                                            '_er-capped_'))
    er.to_zarr(fout, consolidated = True, mode = 'w')

    #Convert sphy and lphy from mol C / m^3 to g C / m^3
    sphy = sphy*12.0107
    lphy = lphy*12.0107

    #Calculate a and b in log B (log10 abundance) vs. log M (log10 gww)
    #in log10 gww
    midsmall = np.log10((mmin+mmid)/2) 
    midlarge = np.log10((mmid+mmax)/2)

    #convert to log10 (gww/size class median size) for log10 abundance
    small = np.log10((sphy*10)/(10**midsmall))
    #convert to log10 (gww/size class median size) for log10 abundance
    large = np.log10((lphy*10)/(10**midlarge))

    #Calculating lope
    slope = (small-large)/(midsmall-midlarge)
    slope.name = 'slope'

    #Creating file path to save slope
    [fout] = [os.path.basename(f) for f in file_list if '_slope_' in f]
    fout = os.path.join(folder_path, '025deg', fout.replace('_slope_', 
                                                            '_slope-capped_'))
    slope.to_zarr(fout, consolidated = True, mode = 'w')

    #intercept is really log10(intercept), same a when small, midsmall are used
    intercept = large-(slope*midlarge)
    intercept.name = 'intercept'

    #Creating file path to save intercept
    [fout] = [os.path.basename(f) for f in file_list if '_intercept_' in f]
    fout = os.path.join(folder_path, '025deg', fout.replace('_intercept_', 
                                                            '_intercept-capped_'))
    intercept.to_zarr(fout, consolidated = True, mode = 'w')

    if experiment in ['obsclim', 'spinup']:
        ui0 = 10**intercept
        ui0.name = 'ui0'
        if experiment == 'obsclim':
            [fout] = glob(os.path.join(folder_path+'_params', '025deg', 
                                       f'ui0_[0-9]*'))
        else:
            [fout] = glob(os.path.join(folder_path+'_params', '025deg',
                                       f'ui0_spinup*'))
        ui0.to_zarr(fout.replace('ui0_', 'ui0-capped_'), consolidated = True, 
                    mode = 'w')


# Merge DBPM outputs to create single yearly or decadal files
# Merging DBPM gridded outputs (yearly or decadal) ----
def merge_files(var, folder, merge_by = 'decade', **kwargs):
    #Getting list of files in folder
    file_list = sorted(glob(os.path.join(folder, f'{var}*')))

    #Getting the base file name to save merged data
    base_file = os.path.basename(file_list[-1])

    #Get years if provided
    if 'years' in kwargs:
        yrs = kwargs.get('years')
    else:
        #Get year from file names
        yrs = np.unique([int(re.findall('_(\\d{4})-\\d{2}', 
                                        os.path.basename(f))[0]) for f in file_list])
    
    #Merge by year or decade
    if merge_by == 'decade':
        dec = np.unique([np.floor(y/10).astype('int') for y in yrs])
        for d in dec:
            #Get files included in current decade
            sub_list = [f for f in file_list if str(d) in f]
            #Getting minimum and maximum years in the decade
            ymax = max([y for y in yrs if str(d) in str(y)])
            ymin = min([y for y in yrs if str(d) in str(y)])
            #Create file name to save merged data
            f_out = re.sub('\\d{4}-\\d{2}|\\d{4}', f'{str(ymin)}-{str(ymax)}', base_file)
            #Open dataset
            ds = xr.open_mfdataset(sub_list, concat_dim = 'time', combine = 'nested')
            #Save dataset
            ds.to_netcdf(os.path.join(folder, f_out))
    if merge_by == 'year':
        for y in yrs:
            #Get files included in current decade
            sub_list = [f for f in file_list if str(y) in f]
            #Create file name to save merged data
            f_out = re.sub('\\d{4}-\\d{2}|\\d{4}', f'{str(y)}', base_file)
            #Open dataset
            ds = xr.open_mfdataset(sub_list, concat_dim = 'time', combine = 'nested')
            #Save dataset
            ds.to_netcdf(os.path.join(folder, f_out))
    else:
        print('"merge_by" should be either "decade" or "year"')

# Defining attributes for coordinates ----
def coord_atts(choice, **extra):
    lat_atts = {'long_name': 'latitude',
                'units': 'degrees_north',
                'standard_name': 'latitude'}
    
    lon_atts = {'long_name': 'longitude',
                'units': 'degrees_east',
                'standard_name': 'longitude'}
    if choice == 'lat':
        lat_atts.update(extra)
        return lat_atts
    elif choice == 'lon':
        lon_atts.update(extra)
        return lon_atts


# Defining attributes for weight bins ----
def weight_atts(choice, **extra):
    # Weight bins
    bin_atts = {'long_name': 'Size-structured bins',
                'units': 'log10 g',
                'standard_name': 'size_structured_weight_bins', 
                'bounds': 'bin_bnds'}
    # Weight bands
    bin_bnds_atts = {'long_name': 
                     'Size-structured bins provide biomass in weight bins: [1-10g],' +
                     ' (10-100g], (100g-1kg], (1-10kg], (10-100kg], (>100kg]',
                     'units': 'g',
                     'standard_name': 'size_structured_weight_bins_range'}
    if choice == 'weights':
        bin_atts.update(extra)
        return bin_atts
    elif choice == 'weight_bands':
        bin_bnds_atts.update(extra)
        return bin_bnds_atts


# Define global attributes ----
def global_atts(**extra):
    global_atts = {'Conventions': 'CF-1.13',
                   'title': 'Dynamic Benthic Pelagic Model (DBPM) outputs: Evaluating ' + 
                   'multiple sources of uncertainty in a global marine ecosystem model ' + 
                   'from a Southern Ocean perspective',
                   'history': f'Created on {dt.now()} by 06_running_gridded_DBPM.py > ' + 
                   '08_calculating_biomass > 14_saving_dbpm_outputs.ipynb',
                   'institution': 'Institute for Marine and Antarctic Studies (IMAS) at ' + 
                   'the University of Tasmania, Australian Centre for Excellence in ' + 
                   'Antarctic Research (ACEAS)',
                   'source': 'DBPM Southern Ocean code repository: ' + 
                   'https://github.com/Benthic-Pelagic-Size-Spectrum-Model/dbpm_southern_ocean',
                   'references': 'Pre-print: http://ssrn.com/abstract=7355654',
                   'comment': 'DBPM simulations followed the FishMIP2.0 ISIMIP3a protocol: ' + 
                   'https://github.com/Fish-MIP/FishMIP2.0_ISIMIP3a',
                   'contact': 'Denisse Fierro Arcos <lilian.fierroarcos@utas.edu.au>, ' + 
                   'Julia Blanchard <julia.blanchard@utas.edu.au>',
                   'license': 'CC-BY-SA-4.0'}
    global_atts.update(extra)
    return global_atts


# Define time attributes ----
def time_atts(**extra):
    time_atts = {'standard_name': 'time',
                 'calendar': 'standard'}
    time_atts.update(extra)
    return time_atts


# Define pelagics attributes ----
def pelagics_atts(choice):
    if choice == 'biomass':
        tpb_atts = {'short_name': 'tpb',
                    'standard_name': 'total_pelagic_biomass_density',
                    'long_name': 
                    'All pelagic consumers (trophic level >1, vertebrates and invertebrates)',
                    'units': 'g m-2',
                    'cell_methods': 'bins: sum'}
        return tpb_atts
    elif choice == 'catches':
        tpc_atts = {'short_name': 'tpc',
                    'standard_name': 'total_pelagic_catch_density_all_fishing_sectors',
                    'long_name': 'Catch at sea of all pelagic consumers (trophic level >1, ' + 
                    'vertebrates and invertebrates)',
                    'units': 'g m-2',
                    'cell_methods': 'bins: sum'}
        return tpc_atts


# Define demersals attributes ----
def demersals_atts(choice):
    if choice == 'biomass':
        tdb_atts = {'short_name': 'tdb',
                    'standard_name': 'total_demersal_biomass_density',
                    'long_name': 
                    'All demersal consumers (trophic level >1, vertebrates and invertebrates)',
                    'units': 'g m-2',
                    'cell_methods': 'bins: sum'}
        return tdb_atts
    elif choice == 'catches':
        tdc_atts = {'short_name': 'tdc',
                    'standard_name': 'total_demersal_catch_density_all_fishing_sectors',
                    'long_name': 'Catch at sea of all demersal consumers (trophic level >1, ' + 
                    'vertebrates and invertebrates)',
                    'units': 'g m-2',
                    'cell_methods': 'bins: sum'}
        return tdc_atts


# Define total community attributes ----
def community_atts(choice):
    if choice == 'biomass':
        tcb_atts = {'short_name': 'tcb',
                    'standard_name': 'total_consumer_biomass_density',
                    'long_name': 
                    'All consumers (trophic level >1, vertebrates and invertebrates)',
                    'units': 'g m-2',
                    'cell_methods': 'bins: sum'}
        return tcb_atts
    elif choice == 'catches':
        tc_atts = {'short_name': 'tc',
                   'standard_name': 'total_catch_density_all_fishing_sectors',
                    'long_name': 'Catch at sea (all catch as a result of all effort including ' + 
                   'reported and IUU) summed for both Industrial and Artisanal sectors',
                    'units': 'g m-2',
                    'cell_methods': 'bins: sum'}
        return tc_atts
        

# Define size-structured community attributes ----
def size_str_comm_atts(choice):
    if choice == 'biomass':
        tcblog10_atts = {'short_name': 'tcblog10',
                         'standard_name': 'total_consumer_biomass_density_in_log10_weight_bins',
                         'long_name': 'Total Consumer Biomass Density in log10 Weight Bins',
                         'units': 'g m-2'}
        return tcblog10_atts
    elif choice == 'catches':
        tclog10_atts = {'short_name': 'tclog10',
                        'standard_name': 
                        'total_catch_density_all_fishing_sectors_in_log10_weight_bins',
                        'long_name': 'Catch at sea (all catch as a result of all effort including ' +
                        'reported and IUU) summed for both Industrial and Artisanal sectors in ' + 
                        'log10 Weight Bins',
                        'units': 'g m-2'}
        return tclog10_atts
        

# Create land mask from outputs
def land_mask(da):
    has_data = da.notnull().any(dim = ['lat', 'lon'])
    first_bin = has_data.argmax('bins').dropna('time')
    mask = (da.sel(time = first_bin.isel(time = 0).time).
        isel(bins = first_bin.isel(time = 0).values))
    mask = xr.where(np.isfinite(mask), 1, np.nan).squeeze().drop_vars(('bins', 'time'))
    return mask


# Applying CF convention to dataset ----
def open_cf_xarray(file_path):
    da = xr.open_mfdataset(file_path)
    da = da.rename_dims({'size_class': 'bins'}).rename_vars({'size_class': 'bins'})
    if da.bins.max() < 100:
        da['bins'] = 10**da.bins
    da['lat'].attrs = coord_atts('lat')
    da['lon'].attrs = coord_atts('lon')
    da['bins'].attrs = weight_atts('weights')
    return da


# Calculating total biomass/catches per functional group (across all weight bins) ----
def total_func_groups(ds_var_int, choice):
    [det_var] = [v for v in ds_var_int.data_vars if 'det' in v]
    [pred_var] = [v for v in ds_var_int.data_vars if 'pred' in v]

    det_atts = demersals_atts(choice)
    pred_atts = pelagics_atts(choice)
    total_atts = community_atts(choice)

    # Create land mask   
    mask = land_mask(ds_var_int[det_var])
    
    #Detritivores
    det = ds_var_int[det_var].sum('bins').to_dataset(name = det_atts['short_name'])
    det.attrs = global_atts()
    det[det_atts['short_name']].attrs = det_atts

    #Predators
    pred = ds_var_int[pred_var].sum('bins').to_dataset(name = pred_atts['short_name'])
    pred.attrs = global_atts()
    pred[pred_atts['short_name']].attrs = pred_atts

    #Community
    total = ((pred[pred_atts['short_name']]+det[det_atts['short_name']]).
        to_dataset(name = total_atts['short_name']))
    total.attrs = global_atts()
    total[total_atts['short_name']].attrs = total_atts
    return det*mask, pred*mask, total*mask


# Calculating total biomass grouped by weight bins ----
def total_func_groups_weight_bins(ds_var_int, weight_range, bins, choice):
    # Create land mask
    mask = land_mask(ds_var_int[list(ds_var_int.keys())[0]])
    
    dict_atts = size_str_comm_atts(choice)
    
    # Format weight bins correctly
    weight_range = weight_range[weight_range.max(axis = 1) > ds_var_int.bins.min().values]
    if len(weight_range)+1 != len(bins):
        bins = bins[bins >= weight_range.min()]
    weight_bins = pd.cut(ds_var_int.bins, bins = bins, labels = weight_range.mean(axis = 1))
    totlog10 = ds_var_int.copy()
    totlog10['bins'] = np.log10(weight_bins.astype('float'))
    totlog10 = (totlog10.groupby('bins').sum().to_array(dim = 'var').sum('var').
        to_dataset(name = dict_atts['short_name']))
    totlog10['bins'].attrs = weight_atts('weights')
    totlog10.attrs = global_atts()
    totlog10['bin_bnds'] = (('bins', 'bnds'), weight_range)
    totlog10['bin_bnds'].attrs = weight_atts('weight_bands')
    totlog10[dict_atts['short_name']].attrs = dict_atts
    return totlog10*mask


# Splitting data over time before saving ----
def split_periods(ds_var_int, period_info, calendar = 'standard'):
    da_split = ds_var_int.sel(time = period_info['dates'])
    da_split['time'] = [cftime.date2num(d, units = period_info['units'], 
                                        calendar = calendar) for d in period_info['dates']]
    da_split['time'].attrs = time_atts()
    da_split['time'] = da_split['time'].assign_attrs(units = period_info['units'])
    return da_split


# Formatting time dimension before saving ----
def format_time(ds_var_int, calendar = 'standard'):
    ds_var_int['time'] = [cftime.date2num(
        d, units = 'days since 1901-1-1 00:00:00', 
        calendar = calendar) for d in pd.DatetimeIndex(ds_var_int.time)]
    ds_var_int['time'].attrs = time_atts()
    ds_var_int['time'] = ds_var_int['time'].assign_attrs(units = 'days since 1901-1-1 00:00:00')
    return ds_var_int
