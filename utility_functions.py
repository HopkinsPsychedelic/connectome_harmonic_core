#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 18:39:04 2020

@author: patricktaylor
"""

import sklearn.neighbors as skn
import time
import numpy as np
import os 
from scipy import sparse

def unmask_medial_wall(masked_feature,medial_wall_mask):
    unmasked_feature=np.zeros(len(medial_wall_mask))
    keepinds=np.where(medial_wall_mask==0)[0]
    unmasked_feature[keepinds]=masked_feature
    
    return unmasked_feature

def unmask_medial_wall_vecs(masked_vecs,medial_wall_mask):
    vecs=np.zeros((len(medial_wall_mask),len(masked_vecs[0,:])))
    for i in range (len(masked_vecs[0,:])):
        vecs[:,i]=unmask_medial_wall(masked_vecs[:,i],medial_wall_mask)
    return vecs

def mask_medial_wall(unmasked_feature,medial_wall_mask):
    keepinds=np.where(medial_wall_mask==0)[0]
    return unmasked_feature[keepinds]

def mask_timeseries(timeseries_unmasked,medial_wall_mask):
    masked_timeseries=np.zeros((len(np.where(medial_wall_mask==0)[0]),len(timeseries_unmasked[0,:])))
    for i in range (len(timeseries_unmasked[0,:])):
        masked_timeseries[:,i]=mask_medial_wall(timeseries_unmasked[:,i],medial_wall_mask)
    return masked_timeseries

def mask_connectivity_matrix(matrix,medial_wall_mask):
    keep_inds=np.where(medial_wall_mask==0)[0]
    return matrix[keep_inds][:,keep_inds]
      


def neighbors(searchedset,queryset,num):
    '''
    computes num nearest neighbors of queryset in searchedset and returns numpy arrays size (len(queryset),num) 
    of indices of searched set and distances between neighbors
    '''
    start=time.time()
    nbrs = skn.NearestNeighbors(n_neighbors=num, algorithm='auto').fit(searchedset)
    distances, indices = nbrs.kneighbors(queryset)
    end=time.time()
    print('[CHAP] Neighbors time =',(end-start))
    return indices,distances

def unmask_medial_wall(masked_feature,medial_wall_mask):
    unmasked_feature=np.zeros(len(medial_wall_mask))
    keepinds=np.where(medial_wall_mask==0)[0]
    unmasked_feature[keepinds]=masked_feature
    
    return unmasked_feature

def unmask_medial_wall_vecs(masked_vecs,medial_wall_mask):
    vecs=np.zeros((len(medial_wall_mask),len(masked_vecs[0,:])))
    for i in range (len(masked_vecs[0,:])):
        vecs[:,i]=unmask_medial_wall(masked_vecs[:,i],medial_wall_mask)
    return vecs
        

def mask_medial_wall(unmasked_feature,medial_wall_mask):
    keepinds=np.where(medial_wall_mask==0)[0]
    return unmasked_feature[keepinds]

def mask_timeseries(timeseries_unmasked,medial_wall_mask):
    masked_timeseries=np.zeros((len(np.where(medial_wall_mask==0)[0]),len(timeseries_unmasked[0,:])))
    for i in range (len(timeseries_unmasked[0,:])):
        masked_timeseries[:,i]=mask_medial_wall(timeseries_unmasked[:,i],medial_wall_mask)
    return masked_timeseries

def demean_timeseries(timeseries):
    means=np.mean(timeseries,axis=1)
    demeaned_timeseries=np.zeros(np.shape(timeseries))
    for i in range (len(timeseries[0,:])):
        demeaned_timeseries[:,i]=timeseries[:,i]-means
    return demeaned_timeseries

def mask_connectivity_matrix(matrix,medial_wall_mask)  :
    keep_inds=np.where(medial_wall_mask==0)[0]
    return matrix[keep_inds][:,keep_inds]
      
def load_all_files_in_directory_sparse(directory,extension='.npz',return_sum=True):
    listofobjects=[]
    for file in os.listdir(directory):
        if file.endswith(extension):
            m=sparse.load_npz(directory+file)
            listofobjects.append(m)
    sum_mat=sparse.csr_matrix(np.shape(m))
    for mat in listofobjects:
        sum_mat+=mat
    return sum_mat


        
    
            
    
    