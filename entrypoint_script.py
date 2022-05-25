#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 11:56:47 2020

@author: bwinston

"""

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

import os
from glob import glob
import input_output as inout
import utility_functions as uts
import argparse
import hcp_preproc_to_chap as hcp_prep
import numpy as np
import bids_to_ch as bids
import email_python


#user inputs cl arguments separated by spaces. args without dashes are required
parser = argparse.ArgumentParser(description='Connectome Harmonic Analysis Pipeline (CHAP)')
parser.add_argument('output_dir', type = str, help = 'CHAP output directory (path)')
parser.add_argument('analysis_level', type = str, help = 'Participant or group mode. Only participant mode supported for now.')
parser.add_argument('--participant_label', nargs='+', help = 'Participant label(s) (not including sub-). If this parameter is not provided all subjects will be analyzed. Multiple participants can be specified with a space separated list')
parser.add_argument('--derivatives_dir', type = str, help = 'Directory containing fMRIprep, Freesurfer, MRtrix3, and Ciftify output directory. Required for CHAP-BIDS pipeline')
parser.add_argument('--hcp_dir', type = str, help = 'HCP (min) preprocessed data directory. First level should be test and retest folders OR if one session just downloads. If test-retest, downloads go in respective session folders. Required for CHAP-HCP pipeline.')
parser.add_argument('--evecs', type = int, help = 'Number of eigenvectors (harmonics) to compute. Default is 100 (minus first trivial harmonic)')
parser.add_argument('--tol', type = int, help = '(Tolerance) search radius of nearest neighbor search for matching endpoints to surface vertices in mm. Default = 1')
parser.add_argument('--sigma', type = int, help = 'Sigma (smoothing)')
parser.add_argument('--epsilon', type = int, help = 'epsilon (smoothing)')
parser.add_argument('--skip_func', type = bool, help = 'Just find structural harmonics, no spectra.')
parser.add_argument('--diff_pipeline', type = str, help = 'Choices: msmt_5tt pipeline or dhollander pipeline based on bids/mrtrix3_connectome. Choose msmt or dholl. Check github sh files for exact commands used.')
parser.add_argument('--streamlines', type = int, help = 'Number of streamlines in MRtrix tckgen')
#parser.add_argument('--mask_med_wall', type = bool, help = 'Mask out medial wall vertices. Default is True.')
parser.add_argument('--binarize', type = bool, help = 'Binarize structural connectivity matrix. Default is True')
parser.add_argument('--recon_only', type = str, help = 'Only compute reconstruction spectrum')
parser.add_argument('--criticality', type = bool, help='compute the criticality of the spectra across subjects')
parser.add_argument('--mem_mb', type=int, help='set maximum memory usage for CHAP')
parser.add_argument('--debug', type = bool, help = 'When True runs without try/except statement (i.e. crashes when a subject fails)')
parser.add_argument('--email_notif', type = str, help = 'Email address where we should send notification that CHAP finished')
args = parser.parse_args() 

##set default arguments if user doesn't supply
if args.mem_mb:
    uts.limit_memory(args.mem_mb)
#read evecs(harmonics) number, set default to 100
if not args.evecs:
    args.evecs = 100
#read tol number, set default to 2
if not args.tol:
    args.tol = 2
#sigma 
if not args.sigma:
    args.sigma = 3
#epsilon
if not args.epsilon:
    args.epsilon = 0.2
#skip func spectra calculation default false
if not args.skip_func:
    args.skip_func = False
#num streamlines default 10 million
if not args.streamlines:
    args.streamlines = '10000000'
else:
    args.streamlines = str(args.streamlines)
#binarize structural connectivity matrix by default
if not args.binarize:
    args.binarize = True
#dholl
if not args.diff_pipeline:
    args.diff_pipeline = 'dholl'
#calculate criticality
if not args.criticality:
    args.criticality = False
#debug mode False default
if not args.debug:
    args.debug = False
    
#create CHAP output directory
inout.if_not_exist_make(f'{args.output_dir}/chap')
#make hcp intermediate dir
if args.hcp_dir: 
    inout.if_not_exist_make(f'{args.output_dir}/chap_work')

#set empty u dict
global u
u = {}

#load mask
u['mask'] = np.load('/home/neuro/repo/hcp_mask.npy')

#find subjects
subs = []
problematic_subs = [] #for later
if args.participant_label: #user input subjects
    subs = [str(sub) for sub in args.participant_label]
elif args.hcp_dir: #get list of hcp subs from data downloaded
    if os.path.exists(f'{args.hcp_dir}/ses-test'): #test-retest data
        sub_list = os.listdir(f'{args.hcp_dir}/ses-test')
    else: 
        sub_list = os.listdir(args.hcp_dir) #one session
    subs = [sub[:6] for sub in sub_list]
    subs = list(dict.fromkeys(subs))    
else: #all subjects from freesurfer output
    subject_dirs = glob(os.path.join(f'{args.derivatives_dir}/freesurfer', "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
print(f'[CHAP] Using sub(s): {subs}')

#main subjects loop
for sub in subs:
    u[f'{sub}_info'] = {}  #create dict in u for each subjs info
    inout.if_not_exist_make(f'{args.output_dir}/chap/sub-{sub}') #subject chap output folder
    #if HCP, run hcp_prep function
    if args.hcp_dir:
        if args.debug == False: #will continue if one subj has error
            try:
                hcp_prep.hcp_chapper(args, sub, u)
            except:
                print(f'Error occurred during {sub}')
                if args.email_notif:
                    email_python.send_email_notif(subject = f'{sub} failed, chap still running')                                                                                                                   
                problematic_subs.append(sub)
        else: #will stop if a subject doesn't work
            hcp_prep.hcp_chapper(args, sub, u)
    #else, run BIDS method
    else:
        if args.debug == False:
            try:
                bids.bids_chapper(u, args, sub)
            except:
                print(f'Error occurred during {sub}')
                problematic_subs.append(sub)
                if args.email_notif:
                    email_python.send_email_notif(subject = f'{sub} failed, chap still running')                                                                                                                   
        else:
            bids.bids_chapper(u, args, sub)
    print(f'[CHAP] Finished {sub}')
print(f'Problematic subs: {problematic_subs}') 
print('''
          )        (                                                   
   (   ( /(  (     )\ )     (                    (         )     (     
   )\  )\()) )\   (()/(     )\         )         )\  (  ( /(  (  )\ )  
 (((_)((_)((((_)(  /(_))  (((_)  (    (    `  ) ((_)))\ )\())))\(()/(  
 )\___ _((_)\ _ )\(_))    )\___  )\   )\  '/(/(  _ /((_|_))//((_)((_)) 
((/ __| || (_)_\(_) _ \  ((/ __|((_)_((_))((_)_\| (_)) | |_(_))  _| |  
 | (__| __ |/ _ \ |  _/   | (__/ _ \ '  \() '_ \) / -_)|  _/ -_) _` |  
  \___|_||_/_/ \_\|_|      \___\___/_|_|_|| .__/|_\___| \__\___\__,_|  
                                          |_|                                                   
Have a pleasant afternoon.                                                                                                                                                                                                 
''')            
if args.email_notif:
    email_python.send_email_notif()                                                                                                                   
                                                                                                                                     
                                                                                                                                      

