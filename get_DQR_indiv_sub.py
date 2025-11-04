#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get DQR for individual subject

Created on Fri Mar 21 10:05:09 2025

@author: smkelley
"""
#%% Imports
import os
import cedalion
import cedalion.nirs
import cedalion.sigproc.quality as quality
import cedalion.xrutils as xrutils


import xarray as xr
import matplotlib.pyplot as p
import cedalion.plots as plots
from cedalion import units
import numpy as np
import pandas as pd
from math import ceil

import sys
sys.path.append('/projectnb/nphfnirs/s/users/shannon/Code/cedalion-pipeline/workflow/scripts/modules')
import module_plot_DQR as dqr
import module_preprocess as preproc

#%%

root_dir = "/projectnb/nphfnirs/s/datasets/Hyperscanning_controls_2425" # CHANGE to your data directory

dyad_num = '11' # CHANGE
subj = 'sub-21' # CHANGE

task_lst = ['JS', 'JR', 'IR', 'PN']#, 'PN'] # CHANGE
runNums = ['05', '04', '02', '02']#, '01']  # CHANGE

ses = 'dyadic' + dyad_num
subj_for_file = subj + '_ses-' + ses


#stim_lst = ['ST', 'DT']  # CHANGE   - your stim names

cfg_dataset = {
    'root_dir' : root_dir,
    'derivatives_subfolder' : '',
}

cfg_prune = {
    'enable' : True,
    'snr_thresh' : 3, # CHANGE
    'sd_thresh' : [1, 45]*units.mm, # CHANGE  # defines the lower and upper bounds for the source-detector separation that we would like to keep
    'amp_thresh' : [1e-3, 1e7], # CHANGE   # define whether a channel's amplitude is within a certain range
    'perc_time_clean_thresh' : 0.6, 
    'sci_threshold' : 0.6,
    'psp_threshold' : 0.1,
    'window_length' : 5 * units.s, 
    'flag_use_sci' : False,   # CHANGE
    'flag_use_psp' : False   # CHANGE
}



cfg_bandpass = { 
    'fmin' : 0.01 * units.Hz, # CHANGE
    'fmax' : 0.5 * units.Hz  # CHANGE
}
        

cfg_preprocess = {
    'median_filt' : 3, # CHANGE  # set to 1 if you don't want to do median filtering
    'prune' : cfg_prune, 
    'freq_filter' : cfg_bandpass,
    #'flag_do_GLM_filter' : True, # CHANGE
}


#%% Load in data

for task, runNum in zip(task_lst, runNums):
    print(f'Getting DQR for subject {subj}, task {task}, run {runNum}')
    subDir = os.path.join(root_dir, f'{subj}', f'ses-{ses}', 'nirs')
    run_nm = f'{subj_for_file}_task-{task}_run-{runNum}'

    stim_lst = [task]

    # update configs for func
    cfg_hrf = {
        'stim_lst' : stim_lst
    }
    cfg_dataset.update({'task' : [task], 'run' : [runNum]})
    
    snirf_path = os.path.join(subDir, run_nm + '_nirs.snirf' )
    
    # check if the snirf file exists
    if not os.path.exists( snirf_path ):
        print( f"Error: File {snirf_path} does not exist" )
    else:
        records = cedalion.io.read_snirf( snirf_path ) 
        rec = records[0]
    
    events_path = os.path.join(subDir, run_nm + '_events.tsv' )
    
    # check if the events.tsv file exists
    if not os.path.exists( events_path ):
        print( f"Error: File {events_path} does not exist" )
    else:
        stim_df = pd.read_csv(events_path, sep='\t' )
        rec.stim = stim_df
    
    #
    # Preprocess 
    #
    
    # Preprocess data with median filt
    #rec = preproc.preprocess( rec, cfg_preprocess['median_filt'] )
    
    # Prune channels
    rec, chs_pruned = preproc.pruneChannels( rec, cfg_preprocess['prune'] )
    pruned_chans = chs_pruned.where(chs_pruned != 0.58, drop=True).channel.values # get array of channels that were pruned
    
    
    # Calculate OD 
    # if flag pruned channels is True, then do rest of preprocessing on pruned amp, if not then do preprocessing on unpruned data
    if cfg_preprocess['prune']['enable']:
        rec["od"] = cedalion.nirs.int2od(rec['amp_pruned'])                
    else:
        rec["od"] = cedalion.nirs.int2od(rec['amp'])
        del rec.timeseries['amp_pruned']   # delete pruned amp from time series
    rec["od_corrected"] = rec["od"]    # need to reassign to new rec_str to work w/ code
    
    # Calculate GVTD on pruned data
    amp_masked = preproc.prune_mask_ts(rec['amp'], pruned_chans)  # use chs_pruned to get gvtd w/out pruned data (could also zscore in gvtd func)
    rec.aux_ts["gvtd"], _ = quality.gvtd(amp_masked) 
    
    
    lambda0 = amp_masked.wavelength[0].wavelength.values
    lambda1 = amp_masked.wavelength[1].wavelength.values
    snr0, _ = quality.snr(amp_masked.sel(wavelength=lambda0), cfg_preprocess['prune']['snr_thresh'])
    snr1, _ = quality.snr(amp_masked.sel(wavelength=lambda1), cfg_preprocess['prune']['snr_thresh'])
    
    # Plot
    dqr.plotDQR( rec, chs_pruned, cfg_preprocess, run_nm, cfg_dataset, cfg_hrf )



#rec, chs_pruned, cfg_preprocess, filenm, cfg_dataset, cfg_hrf

    
# %%
