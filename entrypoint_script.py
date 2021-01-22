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

print('''
 __          __  _                                   _                 _  
 \ \        / / | |                                 | |               | | 
  \ \  /\  / /__| | ___ ___  _ __ ___   ___      ___| |__   __ _ _ __ | | 
   \ \/  \/ / _ \ |/ __/ _ \| '_ ` _ \ / _ \    / __| '_ \ / _` | '_ \| | 
    \  /\  /  __/ | (_| (_) | | | | | |  __/_  | (__| | | | (_| | |_) |_| 
     \/  \/ \___|_|\___\___/|_| |_| |_|\___( )  \___|_| |_|\__,_| .__/(_) 
                                           |/                   | |       
                                                                |_|       
''')

import os, shutil
from glob import glob
import input_output as inout
import argparse
import construct_harmonics as cs
import hcp_preproc_to_chap as hcp_prep
import numpy as np
#user inputs cl arguments separated by spaces. args without dashes are required
#for hcp, hcp_dir is required
#for bids pipeline, qsi_dir, surf_dir, and fs_license_file are required
parser = argparse.ArgumentParser(description='Connectome Harmonic Analysis Pipeline (CHAP)')
parser.add_argument('output_dir', type = str, help = 'CHAP output directory (path)')
parser.add_argument('analysis_level', type = str, help = 'Participant or group mode')
parser.add_argument('--participant_label', type = str, help = 'Participant label(s) (not including sub-). If this parameter is not provided all subjects will be analyzed. Multiple participants can be specified with a list')
parser.add_argument('--qsi_dir', type = str, help = 'qsirecon output directory. Required for BIDS pipeline')
parser.add_argument('--surf_dir', type = str, help = 'BIDS-organized Freesurfer output directory. Required for BIDS pipeline.')
parser.add_argument('--fs_license_file', type = str, help = 'Path to Freesurfer license file (including filename). Required for BIDS pipeline.')
parser.add_argument('--fprep_dir', type = str, help = 'BIDS-organized fMRIprep output dir. Functional images should be in T1w/anat space. Always optional.')
parser.add_argument('--evecs', type = int, help = 'Number of eigenvectors (harmonics) to compute. Default is 100')
parser.add_argument('--nnum', type = int, help = 'Number of nearest neighboring surface vertices to assign to each streamline endpoint. Default = 20' )
parser.add_argument('--hcp_dir', type = str, help = 'HCP (min) preprocessed data directory. First level should be test and retest folders, downloads go in respective session folders. Required for HCP pipeline.')
parser.add_argument('--tol', type = int, help = '(Tolerance) search radius of nearest neighbor search for matching endpoints to surface vertices in mm. Default = 3mm')
parser.add_argument('--mask_med_wall', type = bool, help = 'Mask out medial wall vertices. True or False.H')
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
#create CHAP output directory
inout.if_not_exist_make(f'{args.output_dir}/chap')
#set empty u dict
global u
u = {}
#make hcp intermediate dir, load mask
if args.hcp_dir: 
    inout.if_not_exist_make(f'{args.output_dir}/hcp_preproc')
    u['mask'] = np.load(f'{args.hcp_dir}/hcp_mask.npy')
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
    u[f'{sub}_info'] = {}  #create dict in u for each subjs info
    inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}') #subject chap output folder
    #if HCP, run hcp_prep function
    if args.hcp_dir:
        hcp_prep.hcp_chapper(args, sub, u)
    #else, run BIDS/qsi method
    else:      
        cs.qsi_chap(u, args, sub)
    print(f'[CHAP] Finished {sub}')
print('''
  /$$$$$$  /$$   /$$  /$$$$$$  /$$$$$$$                                                    /$$             /$$                     /$$
 /$$__  $$| $$  | $$ /$$__  $$| $$__  $$                                                  | $$            | $$                    | $$
| $$  \__/| $$  | $$| $$  \ $$| $$  \ $$        /$$$$$$$  /$$$$$$  /$$$$$$/$$$$   /$$$$$$ | $$  /$$$$$$  /$$$$$$    /$$$$$$   /$$$$$$$
| $$      | $$$$$$$$| $$$$$$$$| $$$$$$$/       /$$_____/ /$$__  $$| $$_  $$_  $$ /$$__  $$| $$ /$$__  $$|_  $$_/   /$$__  $$ /$$__  $$
| $$      | $$__  $$| $$__  $$| $$____/       | $$      | $$  \ $$| $$ \ $$ \ $$| $$  \ $$| $$| $$$$$$$$  | $$    | $$$$$$$$| $$  | $$
| $$    $$| $$  | $$| $$  | $$| $$            | $$      | $$  | $$| $$ | $$ | $$| $$  | $$| $$| $$_____/  | $$ /$$| $$_____/| $$  | $$
|  $$$$$$/| $$  | $$| $$  | $$| $$            |  $$$$$$$|  $$$$$$/| $$ | $$ | $$| $$$$$$$/| $$|  $$$$$$$  |  $$$$/|  $$$$$$$|  $$$$$$$
 \______/ |__/  |__/|__/  |__/|__/             \_______/ \______/ |__/ |__/ |__/| $$____/ |__/ \_______/   \___/   \_______/ \_______/
                                                                                | $$                                                  
                                                                                | $$                                                  
                                                                                |__/                                                  

Have a pleasant afternoon.                                                                                                                                                                                                 
''')                                                                                                                                    
                                                                                                                                      
                                                                                                                                      

'''
/Users/bwinston/Documents/fMRI/BIDS/hcp_prep_output participant --hcp_dir /Users/bwinston/Documents/fMRI/BIDS/HCP_Preproc  

'''





