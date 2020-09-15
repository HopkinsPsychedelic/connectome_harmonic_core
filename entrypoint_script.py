#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 11:56:47 2020

@author: bwinston

using https://github.com/BIDS-Apps/example/blob/master/run.py 
https://docs.python.org/3/library/argparse.html
for reference
"""

import argparse
import os
from glob import glob
parser = argparse.ArgumentParser(description='CHAP entrypoint script')
parser.add_argument('input_directory', type = str, help = 'BIDS or HCP input directory')
parser.add_argument('output_directory', type = str, help = 'Output directory')
parser.add_argument('data_type', type = str, help = 'BIDS or HCP')
parser.add_argument('--participant_label', type = str, help = 'Participant label (not including sub-).If this parameter is not provided all subjects should be analyzed. Multiple participants can be specified with a space separated list')

args = parser.parse_args() 

if args.data_type == 'BIDS' or 'bids':     
    subjects_to_analyze = []
    # only for a subset of subjects
    if args.participant_label:
        subjects_to_analyze = args.participant_label
    # for all subjects
    else:
        subject_dirs = glob(os.path.join(args.input_directory, "sub-*"))
        subjects_to_analyze = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]
        
#elif args.data_type == 'HCP' or 'hcp':



#/Users/bwinston/Documents/fMRI/BIDS/1305_dataset/raw_data