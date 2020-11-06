#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 11:56:47 2020

@author: bwinston

see https://github.com/BIDS-Apps/example/blob/master/run.py 
https://docs.python.org/3/library/argparse.html
for reference
"""

import argparse
import os
import shutil
from glob import glob
import numpy as np
import subprocess
import numpy as np
import input_output as inout
import argparse
import decomp as dcp
import utility_functions as uts
import construct_harmonics as cs
#user inputs cl arguments separated by spaces
parser = argparse.ArgumentParser(description='Connectome Harmonic Analysis Pipeline (CHAP)')
parser.add_argument('qsi_dir', type = str, help = 'qsirecon output directory')
parser.add_argument('surf_dir', type = str, help = 'BIDS-organized Freesurfer output directory')
parser.add_argument('output_dir', type = str, help = 'CHAP output directory')
parser.add_argument('analysis_level', type = str, help = 'Participant or group mode')
parser.add_argument('fs_license_file', type = str, help = 'Path to Freesurfer license file (including filename)')
parser.add_argument('--participant_label', type = str, help = 'Participant label(s) (not including sub-). If this parameter is not provided all subjects will be analyzed. Multiple participants can be specified with a space separated list')
parser.add_argument('--fprep_dir', type = str, help = 'BIDS-organized fMRIprep output dir. Functional images should be in T1w/anat space')
parser.add_argument('--parc', type = str, help = "path to parcellation file as vtk with %s for hem")
parser.add_argument('--evecs', type = int, help = 'Number of evecs to compute. Default is 100')
parser.add_argument('--nnum', type = int, help = 'Number of nearest neighboring surface vertices to assign to each streamline endpoint' )
args = parser.parse_args() 
#place Freesurfer license file in freesurfer home dir
#shutil.copyfile(args.fs_license_file, '/opt/freesurfer-6.0.0/license.txt')
#read evecs number, set default to 100
if not args.evecs:
    args.evecs = 100
#read nnum number, set default to 20
if not args.nnum:
    args.nnum = 20
#set empty dict and list
user_info = {}
subs = []
#create output directory
if not os.path.exists(f'{args.output_dir}/chap'): 
    os.mkdir(f'{args.output_dir}/chap') 

#populate dicts with tck output files (reconstruction) for each session, reconstruct surfaces
print('[CHAP] Loading subjects...')
if args.participant_label: #user input subjects
    subs = args.participant_label.split(" ")
else: #for all subjects
    subject_dirs = glob(os.path.join(args.qsi_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]   
for sub in subs:
    user_info[f'{sub}_info'] = {}  #create dict in user_info for each subjs info
    user_info[f'{sub}_info']['streamlines'] = [] #streamlines file locations
    print(f'[CHAP] Reconstructing surfaces for {sub}...')
    os.system(f'bash /home/neuro/repo/surface_resample.sh {args.surf_dir}/sub-{sub}/surf /home/neuro/repo') #run surface reconstruction script
    if not os.path.exists(f'{args.output_dir}/chap/sub-{sub}'): #create output subject folder
        os.mkdir(f'{args.output_dir}/chap/sub-{sub}') 
    if args.fprep_dir:
        user_info[f'{sub}_info']['func'] = [] #functional file locations
    if any('ses' in x for x in os.listdir(f'{args.qsi_dir}/sub-{sub}')): #if multiple sessions
        multises = True
        print(f'[CHAP] Detected multiple sessions for {sub}')
        for ses in os.listdir(f'{args.qsi_dir}/sub-{sub}'): 
            if 'ses' in ses:
                print(f'[CHAP] Locating files for {ses}...')
                os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}') #create output session folders
                for file in os.listdir(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi'):
                    if 'tck' in file:
                        user_info[f'{sub}_info']['streamlines'].append([ses, file]) #streamlines list with each session's .tck
                        print('[CHAP] Located streamlines')
                if args.fprep_dir:
                    print(f'[CHAP] Detected functional images for {sub}')
                    for file in os.listdir(f'{args.fprep_dir}/sub-{sub}/{ses}/func'):
                        if f'space-T1w_desc-preproc_bold.nii.gz' in file:
                            user_info[f'{sub}_info']['func'].append([file])                                    
        for ses, file in user_info[f'{sub}_info']['streamlines']:
            cs.construct_harmonics_calculate_spectra(args, sub, args.output_dir, file, user_info, multises, ses)
    else: #if sub has just one session
        print('[CHAP] Detected only one session')
        for file in os.listdir(f'{args.qsi_dir}/sub-{sub}/dwi'):
            if 'tck' in file:
                user_info[f'{sub}_info']['streamlines'].append([file])
                print('Located streamlines')
        if args.fprep_dir:
            for file in os.listdir(f'{args.fprep_dir}/sub-{sub}/func'):
                if 'space-T1w_desc-preproc_bold.nii.gz' in file: 
                    user_info[f'{sub}_info']['func'].append(file) #functional file locations  
        cs.construct_harmonics_calculate_spectra(sub, args.output_dir, file = user_info[f'{sub}_info']['streamlines'][0])
    


    
  
'''
run config
/Users/bwinston/Documents/fMRI/BIDS/test/qsirecon /Users/bwinston/Documents/fMRI/BIDS/test/freesurfer /Users/bwinston/Documents/fMRI/BIDS/test/output/ participant blah --fprep_dir /Users/bwinston/Documents/fMRI/BIDS/test/fmriprep --participant_label 105923


for hem in ['rh','lh']:
            user_info[f'{sub}_info'][f'{hem}_surf'] = []
            user_info[f'{sub}_info'][f'{hem}_surf'].append(f'{args.surf_dir}/sub-{sub}/surf/{hem}.white') 
            
            
funclist ,img_list = [], []
                    for file in os.listdir(f'{args.fprep_dir}/sub-{sub}/{ses}/func'):
                        if f'space-T1w_desc-preproc_bold.nii.gz' in file:
                            funclist.append(file)
                    for img in funclist:
                        taskstart = img.find('task') + 5
                        img_list.append(img[taskstart:])
                    for fname in img_list:
                        tasklist.append(fname.split('_')[0])
                    tasklist = list(dict.fromkeys(tasklist)) #delete duplicate tasknames
                    for task in tasklist:
                        user_info[f'{sub}_info']['func'][task] = [] #list for each task
                        for img in funclist:
                            if task in img:
                                user_info[f'{sub}_info']['func'][task].append(img)   
'''

















