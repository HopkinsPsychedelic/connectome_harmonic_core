#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 11:56:47 2020

@author: bwinston

see https://github.com/BIDS-Apps/example/blob/master/run.py 
https://docs.python.org/3/library/argparse.html
for reference
test
"""

import argparse
import os
from glob import glob
import numpy as np
import subprocess
import numpy as np
import input_output as inout
import matrix_methods as mm
import argparse
import decomp as dcp
import utility_functions as uts
import compute_spectra as cs
#user inputs cl arguments separated by spaces
parser = argparse.ArgumentParser(description='Connectome Harmonic Analysis Pipeline (CHAP)')
parser.add_argument('qsi_dir', type = str, help = 'qsirecon input directory')
parser.add_argument('surf_dir', type = str, help = 'please input BIDS-organized Freesurfer output dir here')
parser.add_argument('output_dir', type = str, help = 'Output directory')
parser.add_argument('analysis_level', type = str, help = 'participant or group')
parser.add_argument('--participant_label', type = str, help = 'Participant label(s) (not including sub-). If this parameter is not provided all subjects will be analyzed. Multiple participants can be specified with a space separated list')
parser.add_argument('--fprep_dir', type = str, help = 'please input BIDS-organized fMRIprep output dir here. Functional images should be in fsnative space')
parser.add_argument('--parc', type = str, help = "path to parcellation file as vtk with %s for hem")
parser.add_argument('--number', type = str, help = 'number of evecs to compute')
#parser.add_argument('-ps', 'power_spectra', help="option to include or exclude power spectra generation")
#parser.add_argument('-es', 'energy_spectra', help="option to include or exclude energy spectra generation")
#parser.add_argument('-rs', 'reconstruction_spectrum', help="option to include or exclude reconstruction spectrum generation")

args = parser.parse_args() 

#set empty dict and list
user_info = {}
subs = []
if not os.path.exists(f'{args.output_dir}/chap'):
    os.mkdir(f'{args.output_dir}/chap') #create output directory

#populate dicts with tck output files of reconstruction for each session, and {hem}.white
if args.participant_label: #user input subjects
    subs = args.participant_label.split(" ")
else: #for all subjects
    subject_dirs = glob(os.path.join(args.qsi_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]               
for sub in subs:
    user_info[f'{sub}_info'] = {}  #create dict in user_info for each subjs info
    os.mkdir(f'{args.output_dir}/chap/sub-{sub}') #create output subject folder
    #surface files locations
    for hem in ['rh','lh']:
            user_info[f'{sub}_info'][f'{hem}_surf'] = []
            user_info[f'{sub}_info'][f'{hem}_surf'].append(f'{args.surf_dir}/sub-{sub}/surf/{hem}.white') 
    user_info[f'{sub}_info']['streamlines'] = [] #streamlines file locations
    if args.fprep_dir:
        user_info[f'{sub}_info']['func'] = [] #functional file locations
    if 'ses' in os.listdir(f'{args.qsi_dir}/sub-{sub}')[0]: #if multiple sessions
        for ses in os.listdir(f'{args.qsi_dir}/sub-{sub}'): 
            if 'ses' in ses:
                os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}') #create output session folders
                for file in os.listdir(f'{args.qsi_dir}/sub-{sub}/{ses}/dwi'):
                    if 'tck' in file:
                        user_info[f'{sub}_info']['streamlines'].append([ses, file]) #streamlines list with each session's .tck
                if args.fprep_dir:
                    for file in os.listdir(f'{args.fprep_dir}/sub-{sub}/{ses}/func'):
                        for hem in ['L','R']:
                            if f'space-fsnative_hemi-{hem}_bold.func.gii' in file:
                                user_info[f'{sub}_info']['func'].append([ses, file]) #functional file locations               
    else: #if sub has just one session
        #os.mkdir(f'{args.output_dir}/chap/sub-{sub}/ses')
        for file in os.listdir(f'{args.qsi_dir}/sub-{sub}/dwi'):
            if 'tck' in file:
                user_info[f'{sub}_info']['streamlines'].append([file])
        if args.fprep_dir:
            for file in os.listdir(f'{args.fprep_dir}/sub-{sub}/func'):
                        for hem in ['L','R']:
                            if f'space-fsnative_hemi-{hem}_bold.func.gii' in file:
                                user_info[f'{sub}_info']['func'].append(file) #functional file locations
  #NEED TO MAKE ALL OF BELOW WORK FOR PEOPLE WITH JUST ONE SESSION (I.E. NO SES FOLDER)  
    multises = any('ses' in x for x in user_info[f'{sub}_info']['streamlines']) #check whether multiple sessions, set output var
    if multises:
        output = f'{args.output_dir}/chap/sub-{sub}/{ses}'
    else:
        output = f'{args.output_dir}/chap/sub-{sub}'
    for ses, file in user_info[f'{sub}_info']['streamlines']:
        #convert streamlines to .vtk using mrtrix
        tck_name = file.split('/')[-1][:-4]
        os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/endpoints')
        subprocess.check_call("./mrtrix_qsi_pipeline.sh %s %s %s" %(f'{args.qsi_dir}/sub-{sub}/{ses}-test/dwi', tck_name, f'{args.output_dir}/chap/sub-{sub}/{ses}/endpoints'), shell=True)
        #construct surface coordinates, surface endpoints
        lh_surf_path =  user_info[f'{sub}_info']['lh_surf'][0]
        rh_surf_path = user_info[f'{sub}_info']['rh_surf'][0]
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
        #Compute spectra as specified
        #TODO: add correct filepaths once volume-to-surface mapping is complete
        full_path_lh = "placeholder_lh.gii"
        full_path_rh = "placeholder_rh.gii"
        timeseries = cs.read_functional_timeseries(full_path_lh, full_path_rh)
        if(args.fprep_dir):
           os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/powerspectra')
           mean_power_spectrum = cs.mean_power_spectrum(timeseries, vecs)
           dynamic_power_spectrum = cs.dynamic_power_spectrum(timeseries, vecs, vals)
           normalized_power_spectrum = cs.normalized_power_spectrum(timeseries, vecs)
           np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/powerspectra', mean_power_spectrum)
           np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/powerspectra', dynamic_power_spectrum)
           np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/powerspectra', normalized_power_spectrum)
        if(args.fprep_dir):
           os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/energyspectra')
           dynamic_energy_spectrum = cs.dynamic_energy_spectrum(timeseries, vecs, vals)
           normalized_energy_spectrum = cs.normalized_energy_spectrum(timeseries, vecs)
           np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/powerspectra', mean_energy_spectrum)
           np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/powerspectra', dynamic_energy_spectrum)
        if(args.fprep_dir):
            os.mkdir(f'{args.output_dir}/chap/sub-{sub}/{ses}/reconspectra')
            recon_spectrum = cs.dynamic_reconstruction_spectrum(timeseries, vecs, vals)
            np.save(f'{args.output_dir}/chap/sub-{sub}/{ses}/reconspectra', dynamic_reconstruction_spectrum)

        





'''
run config
/Users/bwinston/Documents/fMRI/BIDS/test/qsirecon /Users/bwinston/Documents/fMRI/BIDS/test/freesurfer /Users/bwinston/Documents/fMRI/BIDS/test/output/ participant

'''



hi = ['ricardo', 'winston', 'ses-fasrrre']

multises = any('ses' in x for x in hi)
print(multises)

















