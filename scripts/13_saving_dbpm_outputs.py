#!/usr/bin/env python3

# Loading libraries
import os
from glob import glob
import re
import xarray as xr
import useful_functions as uf
import dask
from datetime import datetime as dt
import pandas as pd
import numpy as np
from itertools import product
from distributed import Client, LocalCluster
from multiprocessing import freeze_support

# This script processes DBPM outputs from the Southern Ocean simulations and
# makes them FishMIP compliant. CF convention is also applied to outputs

# Parent folder where all input variables are stored
base_folder = '/g/data/vf71/la6889/dbpm_inputs/*/gridded_dbpm_outputs/*deg'
base_dirs = glob(base_folder)
runs = ['_simask', '_ccamlr_eff', '']
folder_out = '/g/data/vf71/la6889/dbpm_southern_ocean/outputs_publication/'
combos = list(product(base_dirs, runs, ['biomass', 'catches']))

#Spinup and simulation periods
spinup_period = {'dates': pd.date_range('1841-01-01', '1960-12-31', freq = 'MS'),
                 'units': 'days since 1841-1-1 00:00:00'}
simulation_period = {'dates': pd.date_range('1961-01-01', '2010-12-31', freq = 'MS'),
                     'units': 'days since 1901-1-1 00:00:00'}

# Settting weight bins
bins = np.array([.99, 10, 100, 1e3, 1e4, 1e5, 1e6])
weight_range = np.array([[1, 10], [10, 100], [100, 1000], [1000, 1e4], [1e4, 1e5], [1e5, 1e6]])


def save_outputs(variables, dict_var, period, scenario, region, folder_out, start_sim = '1961'):
    compression_encoding = {'zlib': True, 'complevel': 5, '_FillValue': -9999}
    [var] = [f for f in dict_var[variables].data_vars if 'bin' not in f]
    if period == 'spin':
        fn = f'dbpm_gfdl-mom6-cobalt2_obsclim_histsoc_{scenario}_{var}_{region}_monthly_1841_1960.nc'
    elif period == 'simulation':
        fn = f'dbpm_gfdl-mom6-cobalt2_obsclim_histsoc_{scenario}_{var}_{region}_monthly_{start_sim}_2010.nc'
    # Use a local threaded scheduler inside this worker task, so the dask
    # operations in netcdf_to_zarr do not submit new tasks back to the cluster
    with dask.config.set(scheduler = 'threads'):
        dict_var[variables].to_netcdf(os.path.join(folder_out, fn), format = 'NETCDF4_CLASSIC',
                                      encoding = {var: compression_encoding for var in dict_var[variables].data_vars})


def fishmip_dbpm_outputs(base_dir, run, var_int):
    if '025deg' in base_dir:
        scenario = 'default'
    elif '1deg' in base_dir:
        scenario = '60arcmin'
    reg = re.search('inputs/(.*)/gr', base_dir).group(1)

    fout = os.path.join(folder_out, run.lstrip('_') or 'original')
    os.makedirs(fout, exist_ok = True)
    
    # Get names of relevant files
    if var_int == 'biomass':
        fn = glob(os.path.join(base_dir, f'*_biomass_size_group{run}_fao-*_1841-2010.zarr'))
    elif var_int == 'catches':
        fn = glob(os.path.join(base_dir, f'catches*[t|d]{run}_[0-9]*deg_fao-*_1841*2010*'))
    # Load outputs as CF compliant dataset
    da = uf.open_cf_xarray(fn)
    # Calculate totals for functional groups
    det, pred, total = uf.total_func_groups(da, var_int)
    # Calculate totals for functional groups and size groups
    totlog10 = uf.total_func_groups_weight_bins(da, weight_range, bins, var_int)
    
    if run != '_ccamlr_eff':
        # Create dictionary with spinup results
        spin = {'det': uf.split_periods(det, spinup_period),
                'pred': uf.split_periods(pred, spinup_period),
                'total': uf.split_periods(total, spinup_period),
                'totlog10': uf.split_periods(totlog10, spinup_period)}
        # Create dictionary with simulation results
        simulation = {'det': uf.split_periods(det, simulation_period),
                      'pred': uf.split_periods(pred, simulation_period),
                      'total': uf.split_periods(total, simulation_period),
                      'totlog10': uf.split_periods(totlog10, simulation_period)}
        for v in spin.keys():
            save_outputs(v, spin, 'spin', scenario, reg, fout)
            save_outputs(v, simulation, 'simulation', scenario, reg, fout)
    else:
        # Create dictionary with simulation results
        simulation = {'det': uf.format_time(det),
                      'pred': uf.format_time(pred),
                      'total': uf.format_time(total),
                      'totlog10': uf.format_time(totlog10)}
        for v in simulation.keys():
            save_outputs(v, simulation, 'simulation', scenario, reg, fout, 
                         start_sim = da.time.dt.year.min().values)


#Start cluster
if __name__ == '__main__':
    freeze_support()

    #Start a cluster
    cluster = LocalCluster(threads_per_worker = 4)
    client = Client(cluster)
    client.upload_file(
        '/g/data/vf71/la6889/dbpm_southern_ocean/scripts/useful_functions.py')

    # Process inputs
    futures = client.map(lambda args: fishmip_dbpm_outputs(*args), combos)
   
    # Block until all tasks are done (raises an error if any task fails)
    res = client.gather(futures)
        
    # Shut down cluster cleanly
    client.close()
    cluster.close()
    
