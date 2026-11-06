#!/usr/bin/env python3

# Running gridded DBPM
# Authors: Denisse Fierro Arcos
# Date last update: 2025-08-22

#Loading libraries
import os
from glob import glob
import pandas as pd
import xarray as xr
import json
import numpy as np

# Name of region and model resolution ----
region = 'fao-88'
reg_name = 'west_antarctica'
model_res = '1deg'

# Choose between '_simask' and '_ccamlr_eff' ----
runs = ''#'_ccamlr_eff'

# Paths to input and output folders
base_folder = f'/g/data/vf71/la6889/dbpm_inputs/{reg_name}'
gridded_inputs = os.path.join(base_folder, 'gridded_params', model_res)
if runs == '_simask':
    gridded_outputs = os.path.join(base_folder, 'run_fishing_seaicemask', 
                                   model_res)
elif runs == '_ccamlr_eff':
    gridded_outputs = os.path.join(base_folder, 'run_fishing_new_eff', 
                                   model_res)
else:
    gridded_outputs = os.path.join(base_folder, 'run_fishing', 
                                   model_res)
outputs_folder = os.path.join(base_folder, 'gridded_dbpm_outputs', model_res)
#Ensure outputs folder exists
os.makedirs(outputs_folder, exist_ok = True)

## Load inputs ---
#Gridded DBPM parameters
gridded_params = json.load(open(
    os.path.join(gridded_inputs, f'dbpm_gridded_size_params_{region}_python.json')))

#Mortality from fishing (predators and detritivores)
fish_mort_det = xr.open_zarr(glob(
    os.path.join(gridded_inputs, 'fish-mort-det*'))[0])['fish_mort_det']
fish_mort_pred = xr.open_zarr(glob(
    os.path.join(gridded_inputs, 'fish-mort-pred*'))[0])['fish_mort_pred']

#Size class bins
log10_size_bins_mat = xr.open_zarr('../outputs/log10_size_bins_matrix.zarr/')['size_bins']
size_bin_vals = 10**log10_size_bins_mat

#Area - to be used for masking land areas
area = xr.open_zarr(glob(os.path.join(base_folder, 'gridded', model_res, 
                                      '*areacello*'))[0])['cellareao']

#Effort
effort_init = xr.open_zarr(glob(
    os.path.join(gridded_inputs, 'effort_spinup*'))[0])['effort'].isel(time = 0)

effort = xr.open_mfdataset(glob(
    os.path.join(gridded_outputs, 'effort*[0-9][0-9][0-9].nc')))['effort']

#Add initial effort values
effort = xr.concat([effort_init, effort], dim = 'time')

#Biomass estimates for predators and detritivores
predators = xr.open_mfdataset(glob(
    os.path.join(gridded_outputs, 'predators*[0-9][0-9][0-9].nc')))['predators']

detritivores = xr.open_mfdataset(glob(
    os.path.join(gridded_outputs, 'detritivores*[0-9][0-9][0-9].nc')))['detritivores']

## Calculate fishing mortality - Keep fishing mortality for size classes that are fished
fishing_mort_pred = ((fish_mort_pred*effort).
    isel(size_class = slice(gridded_params['ind_min_fish_pred'], None)))
fishing_mort_det = ((fish_mort_det*effort).
    isel(size_class = slice(gridded_params['ind_min_fish_det'], None)))

## Calculate catches per time step and size class ---
# Predators
catch_pred = fishing_mort_pred*predators*size_bin_vals
catch_pred.name = 'catch_pred'
catch_pred = catch_pred.chunk({'time': 12}).drop_encoding()
catch_pred.to_zarr(os.path.join(outputs_folder,
                                f'catches_pred{runs}_{model_res}_{region}_1841_2010.zarr'), 
                   mode = 'w', consolidated = True)

# Detritivores
catch_det = fishing_mort_det*detritivores*size_bin_vals
catch_det.name = 'catch_det'
catch_det = catch_det.chunk({'time': 12}).drop_encoding()
catch_det.to_zarr(os.path.join(outputs_folder,
                                f'catches_det{runs}_{model_res}_{region}_1841_2010.zarr'), 
                   mode = 'w', consolidated = True)

## Calculate total catches per group and time step ---
# Predators
tot_catch_pred = (catch_pred*gridded_params['log_size_increase']).sum('size_class').drop_encoding()
tot_catch_pred.name = 'tot_catch_pred'
tot_catch_pred.to_zarr(os.path.join(
    outputs_folder, f'tot_catches_pred{runs}_{model_res}_{region}_1841_2010.zarr'), 
                   mode = 'w', consolidated = True)

# Detritivores
tot_catch_det = (catch_det*gridded_params['log_size_increase']).sum('size_class').drop_encoding()
tot_catch_det.name = 'tot_catch_det'
tot_catch_det.to_zarr(os.path.join(
    outputs_folder, f'tot_catches_det{runs}_{model_res}_{region}_1841_2010.zarr'), 
                   mode = 'w', consolidated = True)

## Calculate total catches ---
total_catch = tot_catch_det+tot_catch_pred
total_catch.name = 'total_catches'
total_catch.to_zarr(os.path.join(
    outputs_folder, f'total_catches{runs}_{model_res}_{region}_1841_2010.zarr'), 
                   mode = 'w', consolidated = True)

## Calculate weighted mean for yearly catches
# Creating weights from grid cell area 
weights = area.fillna(0)

# Empty list to store yearly results
catch_weighted_mean = []

#We will group data by year first and then calculate the weighted mean
for yr, d in total_catch.groupby('time.year'):
    d_weighted = d.weighted(weights)
    dw_mean = d_weighted.mean(('time', 'lat', 'lon')).expand_dims({'year': [yr]})
    catch_weighted_mean.append(dw_mean)
#Concatenate into single data array
catch_weighted_mean = xr.concat(catch_weighted_mean, dim = 'year')

# Rename data array before calculations
catch_weighted_mean.name = 'mean_catch'
# Transforming to data frame
catch_weighted_mean = catch_weighted_mean.to_pandas().reset_index()
# Adding units 
catch_weighted_mean['units'] = 't*km-2*year-1'
# Saving results
catch_weighted_mean.to_parquet(os.path.join(
    outputs_folder, f'mean_year_catch_dbpm{runs}_{model_res}_{region}_1841_2010.parquet'))
