#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 11:56:47 2020

@author: bwinston

see https://github.comk/BIDS-Apps/example/blob/master/run.py 
https://docs.python.org/3/library/argparse.html
for reference
"""

#TODO: move entrypoint qsi crap to prep_harmonics_bids function
#TODO: fix sessions stuff for bids method and hcp method

print('[CHAP] Welcome, chap')

import os, shutil
from glob import glob
import input_output as inout
import argparse
import construct_harmonics as cs
import hcp_preproc_to_chap as hcp_prep
#user inputs cl arguments separated by spaces. args without dashes are required
#for hcp, hcp_dir is required
#for bids pipeline, qsi_dir, surf_dir, and fs_license_file are required
parser = argparse.ArgumentParser(description='Connectome Harmonic Analysis Pipeline (CHAP)')
parser.add_argument('output_dir', type = str, help = 'CHAP output directory')
parser.add_argument('analysis_level', type = str, help = 'Participant or group mode')
parser.add_argument('--participant_label', type = str, help = 'Participant label(s) (not including sub-). If this parameter is not provided all subjects will be analyzed. Multiple participants can be specified with a space separated list')
parser.add_argument('--qsi_dir', type = str, help = 'qsirecon output directory. Required for BIDS pipeline')
parser.add_argument('--surf_dir', type = str, help = 'BIDS-organized Freesurfer output directory. Required for BIDS pipeline.')
parser.add_argument('--fs_license_file', type = str, help = 'Path to Freesurfer license file (including filename)')
parser.add_argument('--fprep_dir', type = str, help = 'BIDS-organized fMRIprep output dir. Functional images should be in T1w/anat space. Optional.')
parser.add_argument('--parc', type = str, help = "path to parcellation file as vtk with %s for hem")
parser.add_argument('--evecs', type = int, help = 'Number of eigenvectors to compute. Default is 100')
parser.add_argument('--nnum', type = int, help = 'Number of nearest neighboring surface vertices to assign to each streamline endpoint' )
parser.add_argument('--hcp_dir', type = str, help = 'HCP min. preprocessed data directory. First level should be test and retest folders, downloads go in respective session folders. Required for HCP pipeline.')
parser.add_argument('--tol', type = int, help = '(tolerance) search radius of nearest neighbor search for matching endpoints to surface vertices')
args = parser.parse_args() 
#place Freesurfer license file in freesurfer home dir
if args.fs_license_file:
    shutil.copyfile(args.fs_license_file, '/opt/freesurfer-6.0.0/license.txt')
#read evecs number, set default to 100
if not args.evecs:
    args.evecs = 100
#read tol number, set default to 3
if not args.tol:
    args.tol = 3
#read nnum number, set default to 20
if not args.nnum:
    args.nnum = 20
#make hcp intermediate dir
if args.hcp_dir: 
    inout.if_not_exist_make(f'{args.output_dir}/hcp_preproc')
#create CHAP output directory
inout.if_not_exist_make(f'{args.output_dir}/chap')
#set empty user_info dict
user_info = {}
#find subjects
subs = []
if args.participant_label: #user input subjects
    subs = args.participant_label.split(" ")
elif args.hcp_dir: #get list of hcp subs from data downloaded
    sub_list = os.listdir(f'{args.hcp_dir}/ses-test')
    subs = [sub[:6] for sub in sub_list]
    subs = list(dict.fromkeys(subs))    
else: #all subjects from qsi output
    subject_dirs = glob(os.path.join(args.qsi_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
print(f'[CHAP] Using sub(s): {subs}')
for sub in subs:
    user_info[f'{sub}_info'] = {}  #create dict in user_info for each subjs info
    inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}') #subject chap output folder
    #if HCP, run hcp_prep function
    if args.hcp_dir:
        hcp_prep.hcp_chapper(args, sub, user_info)
    #else, run BIDS/qsi method
    else:      
        cs.qsi_chap(user_info, args, sub)
    print(f'[CHAP] Finished {sub}')
print('[CHAP] CHAP completed. Have a pleasant day.')
 



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















