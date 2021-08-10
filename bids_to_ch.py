#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 21:10:09 2021

@author: bwinston

"""
import input_output as inout
import utility_functions as uts
import subprocess
import os
import compute_spectra as cs
import construct_harmonics as ch

def bids_chapper(u, args, sub): #saves qsiprep tck to sub_info[streamlines]; passes off to ciftify_chap
    u[f'{sub}_info']['streamlines'] = [] #where streamlines files will go
    if any('ses' in x for x in os.listdir(f'{args.qsi_dir}/sub-{sub}')): #if multiple sessions
        multises = True
        print(f'[CHAP] Detected multiple sessions for {sub}')
        for ses in os.listdir(f'{args.mrtrix_dir}/sub-{sub}'): 
            if 'ses' in ses:
                u[f'{sub}_info'][ses] = {}
                inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}') #create output session folders
                
                
                
                for file in os.listdir(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi'): #look in qsirecon output dir for tck
                    if 'tck' in file:
                        u[f'{sub}_info'][ses]['streamlines'] = file #streamlines list with each session's track file
                        print(f'[CHAP] Located tractography for sub-{sub} {ses}')
            get_endpoints(args, sub, u, multises, ses) 
    else: #if sub has just one session
        print(f'[CHAP] Detected only one session for {sub}')
        multises = False
        ses = '' 
        u[f'{sub}_info'][ses] = {}
        for file in os.listdir(f'{args.qsi_dir}/sub-{sub}/dwi'):
            if 'tck' in file:
                u[f'{sub}_info'][ses]['streamlines'] = file
                print(f'[CHAP] Located tractography for sub-{sub}')
        get_endpoints(args, sub, u, multises, ses)

def get_endpoints(args, sub, u, multises, ses):
    tck_name = u[f'{sub}_info'][ses]['streamlines'].split('/')[-1][:-4]
    if not os.path.exists(f'{args.output_dir}/chap/sub-{sub}/{ses}/{tck_name}_endpoints.vtk'):
        subprocess.check_call("/home/neuro/repo/mrtrix_qsi_pipeline.sh %s %s %s" %(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi', tck_name, f'{args.output_dir}/chap/sub-{sub}/{ses}'), shell=True) #run mrtrix bash script
    u[f'{sub}_info'][ses]['endpoints'] = f'{args.output_dir}/chap/sub-{sub}/{ses}/{tck_name}_endpoints.vtk' #save endpoints to u
    for file in os.listdir(f'{args.output_dir}/chap/sub-{sub}/{ses}'):
        if '_endpoints.tck' in file:
            os.remove(f'{args.output_dir}/chap/sub-{sub}/{ses}/{file}') #remove endpoints tck
    ciftify_chap(u, args, sub, multises, ses) 
        
def ciftify_chap(u, args, sub, multises, ses):
    #get ciftify surfs
    u[f'{sub}_info'][ses]['surfs'] = {}
    u[f'{sub}_info'][ses]['surfs']['lh'] = f'{args.ciftify_dir}/sub-{sub}/T1w/fsaverage_LR32k/sub-{sub}.L.white.32k_fs_LR.surf.gii'
    u[f'{sub}_info'][ses]['surfs']['rh'] = f'{args.ciftify_dir}/sub-{sub}/T1w/fsaverage_LR32k/sub-{sub}.R.white.32k_fs_LR.surf.gii'
    if multises == True:
        print('[CHAP] PLEASE NOTE: You have input a dataset with multiple sessions. Ciftify only calculates one surface, which will be used at multiple sessions. This is not a problem (see Winston et. al 2021)')
    print(f'[CHAP] Found ciftify surfaces for {sub}')
    #for each ses, {sub}_info[ses][func] is a list of the dtseries files
    if os.path.exists(f'{args.ciftify_dir}/sub-{sub}/MNINonLinear/Results'): #there is functional stuff
        u[f'{sub}_info'][ses]['is_func'] = 'cift'
        u[f'{sub}_info'][ses]['func'] = []
        for func_dir in os.listdir(f'{args.ciftify_dir}/sub-{sub}/MNINonLinear/Results'): 
            if ses in func_dir: #ses can be empty, remember
                for file in os.listdir(f'{args.ciftify_dir}/sub-{sub}/MNINonLinear/Results/{func_dir}'): #ciftify functional timeseries directories
                    if 'dtseries' in file:
                        u[f'{sub}_info'][ses]['func'].append(f'{args.ciftify_dir}/sub-{sub}/MNINonLinear/Results/{func_dir}/{file}')
                        print(f'[CHAP] Found ciftify timeseries: {file}') 
    ch.construct_harmonics(args, sub, ses, u, multises) 

def bids_spectra_prep(args,sub,ses,u,vecs,vals):
    for dts in u[f'{sub}_info'][ses]['func']: #each ciftify dtseries
        bids_stuff = f'sub-{sub}_{inout.get_bids_stuff(dts)}' #e.g. sub-{sub}_ses-{ses}_task-{task}
        inout.dts_to_func_gii(dts, f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{bids_stuff}') #extract cortical timeseries with connectome workbench
        u[f'{sub}_info'][ses][f'{bids_stuff}_ts'] = cs.read_functional_timeseries(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{bids_stuff}_hem-l.func.gii', f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{bids_stuff}_hem-r.func.gii') #func.gii to timeseries
        u[f'{sub}_info'][ses][f'{bids_stuff}_ts'] = uts.mask_timeseries(u[f'{sub}_info'][ses][f'{bids_stuff}_ts'], u['mask']) #mask timeseries
        ch.func_spectra(args, sub, ses, u[f'{sub}_info'][ses][f'{bids_stuff}_ts'], inout.get_task(dts), bids_stuff, vecs, vals)
   
    










