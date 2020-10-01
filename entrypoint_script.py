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
import subprocess
import numpy as np
#import input_output as inout
import matrix_methods as mm
import argparse
import decomp as dcp
import utility_functions as uts
#user inputs cl arguments separated by spaces
parser = argparse.ArgumentParser(description='CHAP entrypoint script')
parser.add_argument('--input_dir', type = str, help = 'qsirecon input directory')
parser.add_argument('--surf_dir', type = str, help = 'please input BIDS-organized Freesurfer output dir here')
parser.add_argument('--output_dir', type = str, help = 'Output directory')
parser.add_argument('--analysis_level', type = str, help = 'participant or group')
parser.add_argument('--participant_label', type = str, help = 'Participant label(s) (not including sub-). If this parameter is not provided all subjects will be analyzed. Multiple participants can be specified with a space separated list')
parser.add_argument('-p','--parc',help="path to parcellation file as vtk with %s for hem")
parser.add_argument('-n', '--number', help='number of evecs to compute')
args = parser.parse_args() 

#set empty dict and list
user_info = {}
subs = []
if not os.path.exists(f'{args.output_dir}/chap'):
    os.mkdir(f'{args.output_dir}/chap') #create output directory

#populate dicts with tck output files of reconstruction for each session, and {hem}.white
#user input subjects
if args.participant_label:
    subs = args.participant_label.split(" ")
#for all subjects
else:
    subject_dirs = glob(os.path.join(args.input_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]                 
for sub in subs:
    user_info[f'{sub}_info'] = {}  #create dict in user_info for each subjs info
    os.mkdir(f'{args.output_dir}/chap/sub-{sub}') #create output subject folder
    #surface files locations
    for hem in ['rh','lh']:
            user_info[f'{sub}_info'][f'{hem}_surf'] = []
            user_info[f'{sub}_info'][f'{hem}_surf'].append(f'{args.surf_dir}_sub-{sub}/surf/{hem}.white')
    #streamlines file locations
    user_info[f'{sub}_info']['streamlines'] = []
    if 'ses' in os.listdir(f'{args.input_dir}/sub-{sub}')[0]: #if multiple sessions
        for ses in os.listdir(f'{args.input_dir}/sub-{sub}'): 
            if 'ses' in ses:
                os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}') #create output session folders
                for file in os.listdir(f'{args.input_dir}/sub-{sub}/{ses}-test/dwi'):
                    if 'tck' in file:
                        user_info[f'{sub}_info']['streamlines'].append([ses, file]) #streamlines list with each session's .tck
    else: #if sub has just one session
        os.mkdir(f'{args.output_dir}/chap/sub-{sub}/ses')
        for file in os.listdir(f'{args.input_dir}/sub-{sub}/ses-test/dwi'):
            if 'tck' in file:
                user_info[f'{sub}_info']['streamlines'].append(['ses', file])


    for ses, file in user_info[f'{sub}_info']['streamlines']:
    #convert streamlines to .vtk using mrtrix
        tck_name = file.split('/')[-1][:-4]
        os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/endpoints')
        subprocess.check_call("./mrtrix_qsi_pipeline.sh %s %s %s" %(f'{args.input_dir}/sub-{sub}/{ses}-test/dwi', tck_name, f'{args.output_dir}/chap/sub-{sub}/{ses}/endpoints'), shell=True)
	#construct surface coordinates, surface endpoints
        lh_surf_path =  user_info[f'{sub}_info']['lh_surf']
        rh_surf_path = user_info[f'{sub}_info']['rh_surf']
        if lh_surf_path.endswith('.vtk'):
            sc,si=inout.read_vtk_surface_both_hem(lh_surf_path, rh_surf_path)
        else:
            sc,si=inout.read_gifti_surface_both_hem(lh_surf_path, rh_surf_path)
        streamline_path = f'{args.output_dir}/chap/sub-{sub}/{ses}/endpoints/{tck_name}.vtk'
        ec=inout.read_streamline_endpoints(streamline_path)
        surf_mat=mm.construct_surface_matrix(sc,si)
        ihc_mat=mm.construct_inter_hemi_matrix(sc,tol=3)
        struc_conn_mat=mm.construct_structural_connectivity_matrix(sc,ec,tol=3,NNnum=45)
        vals,vecs=dcp.lapDecomp(struc_conn_mat,args.number)
        os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/vecsvals')
        os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis')
        np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/vecsvals/',[vals,vecs])
        inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/',sc,si,vecs)





'''
run config
/Users/bwinston/Documents/fMRI/BIDS/test/qsirecon /Users/bwinston/Documents/fMRI/BIDS/test/freesurfer /Users/bwinston/Documents/fMRI/BIDS/test/output/ participant

'''




















