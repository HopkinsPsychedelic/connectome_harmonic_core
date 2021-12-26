#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 21:10:09 2021

@author: bwinston

"""
import input_output as inout
import utility_functions as uts
import os
import compute_spectra as cs
import construct_harmonics as ch
import hcp_preproc_to_chap as hcp_prep

def bids_chapper(u, args, sub): #saves qsiprep tck to sub_info[streamlines]; passes off to ciftify_chap
    u[f'{sub}_info']['streamlines'] = [] #where streamlines files will go
    freesurfer_dir = f'{args.freesurfer_dir}/sub-{sub}'
    if any('ses' in x for x in os.listdir(f'{args.mrtrix_dir}/sub-{sub}')): #if multiple sessions
        multises = True
        print(f'[CHAP] Detected multiple sessions for {sub}')
        for ses in os.listdir(f'{args.mrtrix_dir}/sub-{sub}'): 
            if 'ses' in ses:
                u[f'{sub}_info'][ses] = {}
                inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}') #create output session folders               
                diffusion_dir = f'{args.mrtrix_dir}/sub-{sub}/{ses}/dwi/sub-{sub}_ses-{ses}'
                hcp_prep.mrtrix_recon(u,sub, ses, args, f'{diffusion_dir}_desc-preproc_dwi.nii.gz', f'{diffusion_dir}_desc-preproc_dwi.bval', f'{diffusion_dir}_desc-preproc_dwi.bvec', freesurfer_dir, f'{diffusion_dir}_desc-brain_mask.nii.gz')            
            ciftify_chap(u, args, sub, multises, ses) 
    else: #if sub has just one session
        print(f'[CHAP] Detected only one session for {sub}')
        multises = False
        ses = '' 
        u[f'{sub}_info'][ses] = {}
        diffusion_dir = f'{args.mrtrix_dir}/sub-{sub}/{ses}/dwi/sub-{sub}'
        hcp_prep.mrtrix_recon(u,sub, ses, args, f'{diffusion_dir}_desc-preproc_dwi.nii.gz', f'{diffusion_dir}_desc-preproc_dwi.bval', f'{diffusion_dir}_desc-preproc_dwi.bvec', freesurfer_dir, f'{diffusion_dir}_desc-brain_mask.nii.gz')            
        ciftify_chap(u, args, sub, multises, ses)
        
def ciftify_chap(u, args, sub, multises, ses):
    #get ciftify surfs
    u[f'{sub}_info'][ses]['surfs'] = {}
    u[f'{sub}_info'][ses]['surfs']['lh'] = f'{args.ciftify_dir}/sub-{sub}/T1w/fsaverage_LR32k/sub-{sub}.L.white.32k_fs_LR.surf.gii'
    u[f'{sub}_info'][ses]['surfs']['rh'] = f'{args.ciftify_dir}/sub-{sub}/T1w/fsaverage_LR32k/sub-{sub}.R.white.32k_fs_LR.surf.gii'
    if multises == True:
        print('[CHAP] PLEASE NOTE: You have input a dataset with multiple sessions. Ciftify only calculates one surface, which will be used at multiple sessions. This is not a problem (see Winston et. al 2022)')
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
    if os.path.exists(f'{args.output_dir}/chap/sub-{sub}/{ses}/vecs.npy'):
        print('[CHAP] Harmonics already detected. Checking for spectra...')
        check_func(args,sub,ses,u,vecs,vals)
    else:
        ch.construct_harmonics(args, sub, ses, u, multises) 

def bids_spectra_prep(args,sub,ses,u,vecs,vals):
    if 'HCP_Raw' in args.ciftify_dir: #if inputting HCP Raw data to BIDS version (used in Winston et. al 2022)
        tasks = ['WM','MOTOR','LANGUAGE','EMOTION','GAMBLING','SOCIAL','RELATIONAL']
        for task in tasks:
            for dts in u[f'{sub}_info'][ses]['func']: #each ciftify dtseries
                if task in dts:
                    for dire in ['LR','RL']:
                        bids_stuff = f'sub-{sub}_{ses}_task-{task}_acq-{dire}'
                        u[f'{sub}_info'][ses][f'{task}_{dire}'] = dts
                        inout.dts_to_func_gii(u[f'{sub}_info'][ses][f'{task}_{dire}'], f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{bids_stuff}')
                        u[f'{sub}_info'][ses][f'{task}_{dire}'] = cs.read_functional_timeseries(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{bids_stuff}_hem-l.func.gii', f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{bids_stuff}_hem-r.func.gii')
                        u[f'{sub}_info'][ses][f'{task}_{dire}'] = uts.mask_timeseries(u[f'{sub}_info'][ses][f'{task}_{dire}'],u['mask'])
                        os.remove(f'{func_dir}/{bids_stuff}_hem-l.func.gii')
                        os.remove(f'{func_dir}/{bids_stuff}_hem-r.func.gii')
                    print(f'[CHAP] Concatenating LR and RL PE direction scans for {sub} {ses} {task} scan...')
                    u[f'{sub}_info'][ses][f'{task}_ts'] = inout.combine_pe(u[f'{sub}_info'][ses][f'{task}_{LR}'],u[f'{sub}_info'][ses][f'{task}_RL'])  
                    ch.func_spectra(args,sub,ses,u[f'{sub}_info'][ses][f'{task}_ts'],task,bids_stuff,vecs,vals)   
    else: #non HCP raw data
        for dts in u[f'{sub}_info'][ses]['func']: #each ciftify dtseries
            bids_stuff = f'sub-{sub}_{inout.get_bids_stuff(dts)}' #e.g. sub-{sub}_ses-{ses}_task-{task}
            inout.dts_to_func_gii(dts, f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{bids_stuff}') #extract cortical timeseries with connectome workbench
            u[f'{sub}_info'][ses][f'{bids_stuff}_ts'] = cs.read_functional_timeseries(f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{bids_stuff}_hem-l.func.gii', f'{args.output_dir}/chap/sub-{sub}/{ses}/func/{bids_stuff}_hem-r.func.gii') #func.gii to timeseries
            u[f'{sub}_info'][ses][f'{bids_stuff}_ts'] = uts.mask_timeseries(u[f'{sub}_info'][ses][f'{bids_stuff}_ts'], u['mask']) #mask timeseries
            os.remove(f'{func_dir}/{bids_stuff}_hem-l.func.gii')
            os.remove(f'{func_dir}/{bids_stuff}_hem-r.func.gii')
            ch.func_spectra(args, sub, ses, u[f'{sub}_info'][ses][f'{bids_stuff}_ts'], inout.get_task(dts), bids_stuff, vecs, vals)
            
    










