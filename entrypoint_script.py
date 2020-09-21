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
from glob import glob

#user inputs cl arguments separated by spaces
parser = argparse.ArgumentParser(description='CHAP entrypoint script')
parser.add_argument('input_dir', type = str, help = 'qsirecon or HCP input directory')
parser.add_argument('output_dir', type = str, help = 'Output directory')
parser.add_argument('data_type', type = str, help = 'BIDS or HCP')
parser.add_argument('--participant_label', type = str, help = 'Participant label(s) (not including sub-). If this parameter is not provided all subjects will be analyzed. Multiple participants can be specified with a space separated list')
parser.add_argument('--surf_dir', type=str, help = 'If input_dir = qsirecon, please input BIDS-organized Freesurfer output dir here')
args = parser.parse_args() 

#set empty dict and list
user_info = {}
subs = []

#for HCP populate user_info with the files we want (i.e. the two surfaces and whatever diffusion shit that patrick uses)
#need to see how HCP data are organized on JHPCE
if args.data_type == 'HCP':
    #populate subs list
    # only for a subset of subjects
    if args.participant_label:
        subs = args.participant_label
    # for all subjects. Assuming HCP input directory has a bunch of subject folders (i.e. /105923, /113213 etc.)
    else:
        print()
#if data_type is BIDS, takes in qsirecon, freesurfer output dirs, populate dicts with tck output files of reconstruction for each session, and {hem}.white
elif args.data_type == 'BIDS':
    #user input subjects
    if args.participant_label:
        subs = args.participant_label
    #for all subjects
    else:
        subject_dirs = glob(os.path.join(args.input_dir, "sub-*"))
        subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]                 
    for sub in subs:
        user_info[f'{sub}_info'] = {}  #create dict in user_info for each subjs info
        #surface files locations
        for hem in ['rh','lh']:
                user_info[f'{sub}_info'][f'{hem}_surf'] = []
                user_info[f'{sub}_info'][f'{hem}_surf'].append(f'{args.surf_dir}_sub-{sub}/surf/{hem}.white')
        #check if sub has multiple sessions
        if 'ses' in os.listdir(f'{args.input_dir}/sub-{sub}')[0]: #if multiple sessions
            for ses in os.listdir(f'{args.input_dir}/sub-{sub}'): 
                if 'ses' in ses:
                    user_info[f'{sub}_info'][ses] = {} #create a dict for each session
                    #streamlines file location
                    user_info[f'{sub}_info'][ses]['streamlines'] = []
                    for file in os.listdir(f'{args.input_dir}/sub-{sub}/{ses}/dwi'):
                        if 'tck' in file:
                            user_info[f'{sub}_info'][ses]['streamlines'].append(file)
        else: #if sub has just one session
            #streamlines file location
            user_info[f'{sub}_info'][ses]['streamlines'] = [] 
            for file in os.listdir(f'{args.input_dir}/sub-{sub}/dwi'):
                if 'tck' in file:
                    user_info[f'{sub}_info'][ses]['streamlines'].append(file)
            

































