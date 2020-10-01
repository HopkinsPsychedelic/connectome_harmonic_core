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
parser.add_argument('input_dir', type = str, help = 'qsirecon input directory')
parser.add_argument('surf_dir', type = str, help = 'please input BIDS-organized Freesurfer output dir here')
parser.add_argument('output_dir', type = str, help = 'Output directory')
parser.add_argument('analysis_level', type = str, help = 'participant or group')
parser.add_argument('--participant_label', type = str, help = 'Participant label(s) (not including sub-). If this parameter is not provided all subjects will be analyzed. Multiple participants can be specified with a space separated list')
parser.add_argument('-v','--savevis',help='save evecs to surface, if used, give path to output vtk')
parser.add_argument('-p','--parc',help="path to parcellation file as vtk with %s for hem")
parser.add_argument('-n', '--number', help='number of evecs to compute')
parser.add_argument('-o','--outputvecs',help='path to output vector/ vals file')
args = parser.parse_args() 

#set empty dict and list
user_info = {}
subs = []
if not os.path.exists(f'{args.output_dir}/chap'):
    os.mkdir(f'{args.output_dir}/chap') #create output directory

#populate dicts with tck output files of reconstruction for each session, and {hem}.white
#user input subjects
if args.participant_label:
    subs = args.participant_label
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
                for file in os.listdir(f'{args.input_dir}/sub-{sub}/{ses}/dwi'):
                    if 'tck' in file:
                        user_info[f'{sub}_info']['streamlines'].append(file) #streamlines list with each session's .tck
    else: #if sub has just one session
        for file in os.listdir(f'{args.input_dir}/sub-{sub}/dwi'):
            if 'tck' in file:
                user_info[f'{sub}_info']['streamlines'].append(file)




                
'''                
	#convert streamlines to .vtk using mrtrix
    tck_name = user_info[f'{sub}_info']['streamlines'].split('/')[-1][:-5]
    subprocess.check_call("./mrtrix_qsi_pipeline.sh %s %s" %(f'{args.input_dir}/sub-{sub}/dwi', tck_name), shell=True)            
    streamline_path = f'/home/neuro/output/endpoints/{tck_name}.tck'
	#construct surface coordinates, surface endpoints
    lh_surf_path =  user_info[f'{sub}_info']['lh_surf']
    rh_surf_path = user_info[f'{sub}_info']['rh_surf']
    if lh_surf_path.endswith('.vtk'):
            sc,si=inout.read_vtk_surface_both_hem(lh_surf_path, rh_surf_path)
    else:
            sc,si=inout.read_gifti_surface_both_hem(lh_surf_path, rh_surf_path)
    ec=inout.read_streamline_endpoints(streamline_path)
    surf_mat=mm.construct_surface_matrix(sc,si)
    ihc_mat=mm.construct_inter_hemi_matrix(sc,tol=3)
    struc_conn_mat=mm.construct_structural_connectivity_matrix(sc,ec,tol=3,NNnum=45)
    mask=inout.generate_mask_from_parc(args.parc % 'lh', args.parc % 'rh')
    masked_mat = uts.mask_connectivity_matrix(surf_mat+ihc_mat+struc_conn_mat,mask)
    vals,vecs=dcp.lapDecomp(masked_mat,args.number)
    np.save(args.outputvecs,[vals,vecs])
    if args.savevis:
            inout.save_eigenvector(args.savevis,sc,si,vecs)



'''


'''
run config
/Users/bwinston/Documents/fMRI/BIDS/test/qsirecon /Users/bwinston/Documents/fMRI/BIDS/test/freesurfer /Users/bwinston/Documents/fMRI/BIDS/test/output/ participant

'''




















