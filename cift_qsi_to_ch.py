#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 21:10:09 2021

@author: bwinston
"""
import decomp as dcp
import input_output as inout
import utility_functions as uts
import subprocess
import numpy as np
import os
import matrix_methods as mm
import compute_spectra as cs
from scipy import sparse
from itertools import product
import construct_harmonics as ch

def qsi_chap(u, args, sub): #saves qsiprep tck to sub_info[streamlines]; passes off to ciftify_chap
    u[f'{sub}_info']['streamlines'] = [] #where streamlines files will go
    if any('ses' in x for x in os.listdir(f'{args.qsi_dir}/sub-{sub}')): #if multiple sessions
        multises = True
        print(f'[CHAP] Detected multiple sessions for {sub}')
        for ses in os.listdir(f'{args.qsi_dir}/sub-{sub}'): 
            if 'ses' in ses:
                print(f'[CHAP] Locating files for {ses}...')
                u[f'{sub}_info'][ses] = {}
                inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}/{ses}') #create output session folders
                for file in os.listdir(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi'): #look in qsirecon output dir for tck
                    if 'tck' in file:
                        u[f'{sub}_info'][ses]['streamlines'] = file #streamlines list with each session's .tck
                        print('[CHAP] Located streamlines')
            get_endpoints(args, sub, u, multises, ses) 
    else: #if sub has just one session
        print('[CHAP] Detected only one session')
        multises = False
        for file in os.listdir(f'{args.qsi_dir}/sub-{sub}/dwi'):
            if 'tck' in file:
                u[f'{sub}_info']['streamlines'] = file
                print('[CHAP] Located streamlines')
        get_endpoints(args, sub, u, multises, ses = '')

def get_endpoints(args, sub, u, multises, ses):
    tck_name = u[f'{sub}_info'][ses]['streamlines'].split('/')[-1][:-4]
    print('[CHAP] Saving streamline endpoints and converting to vtk...')
    subprocess.check_call("/home/neuro/repo/mrtrix_qsi_pipeline.sh %s %s %s" %(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi', tck_name, f'{args.output_dir}/chap/sub-{sub}/{ses}'), shell=True) #run mrtrix bash script
    u[f'{sub}_info'][ses]['endpoints'] = f'{args.output_dir}/chap/sub-{sub}/{ses}/{tck_name}_endpoints.vtk' #save endpoints to u
    for file in os.listdir(f'{args.output_dir}/chap/sub-{sub}/{ses}'):
        if '_endpoints.tck' in file:
            os.remove(f'{args.output_dir}/chap/sub-{sub}/{ses}/{file}') #remove endpoints tck
    print('[CHAP] Finished MRtrix commands')
    ciftify_chap(u, args, sub, multises) 
        
def ciftify_chap(u, args, sub, multises, ses):
    #get ciftify surfs
    u[f'{sub}_info'][ses]['surfs']['lh'] = f'{args.ciftify_dir}/sub-{sub}/T1w/fsaverage_LR32k/sub-{sub}.L.white.32k_fs_LR.surf.gii'
    u[f'{sub}_info'][ses]['surfs']['rh'] = f'{args.ciftify_dir}/sub-{sub}/T1w/fsaverage_LR32k/sub-{sub}.R.white.32k_fs_LR.surf.gii'
    if multises == True:
        print('[CHAP] PLEASE NOTE: You have input a dataset with multiple sessions. Ciftify only calculates one surface, which will be used at multiple sessions. This is not a problem (see Winston et. al 2022)')
    print(f'[CHAP] Found ciftify surfaces for {sub}')
    #for each ses, sub_info[ses][func] is a list of the dtseries files
    if os.path.exists(f'{args.ciftify_dir}/sub-{sub}/MNINonLinear/Results'):
        u[f'{sub}_info'][ses]['is_func'] == True
        u[f'{sub}_info'][ses]['func'] = []
        for func_dir in os.listdir(f'{args.ciftify_dir}/sub-{sub}/MNINonLinear/Results'):
            if ses in func_dir:
                for file in os.listdir(f'{args.ciftify_dir}/sub-{sub}/MNINonLinear/Results/{func_dir}'):
                    if 'dtseries' in file:
                        u[f'{sub}_info'][ses]['func'].append(file)
                        print(f'[CHAP] Found ciftify timeseries: {file}')
    else:
       u[f'{sub}_info'][ses]['is_func'] == False 
    ch.construct_harmonics_calculate_spectra(args, sub, ses, u, multises) 













