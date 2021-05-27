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
import icc

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

def procrustes(source, target, center=False, scale=False):
    """Align `source` to `target` using procrustes analysis.
    Parameters
    ----------
    source : 2D ndarray, shape = (n_samples, n_feat)
        Source dataset.
    target : 2D ndarray, shape = (n_samples, n_feat)
        Target dataset.
    center : bool, optional
        Center data before alignment. Default is False.
    scale : bool, optional
        Remove scale before alignment. Default is False.
    Returns
    -------
    aligned : 2D ndarray, shape = (n_samples, n_feat)
        Source dataset aligned to target dataset.
    """

    # Translate to origin
    if center:
        ms = source.mean(axis=0)
        mt = target.mean(axis=0)

        source = source - ms
        target = target - mt

    # Remove scale
    if scale:
        ns = np.linalg.norm(source)
        nt = np.linalg.norm(target)
        source /= ns
        target /= nt

    # orthogonal transformation: rotation + reflection
    u, w, vt = np.linalg.svd(target.T.dot(source).T)
    
    
    t = u.dot(vt)
    
    #print(np.shape(t))
    # Recover target scale
    if scale:
        t *= w.sum() * nt

    aligned = source.dot(t)
    if center:
        aligned += mt
    return aligned

def procrustes_alignment(data, reference=None, n_iter=10, tol=1e-5,
                         return_reference=False, verbose=False):
    """Iterative alignment using generalized procrustes analysis.
    Parameters
    ----------
    data :  list of ndarray, shape = (n_samples, n_feat)
        List of datasets to align.
    reference : ndarray, shape = (n_samples, n_feat), optional
        Dataset to use as reference in the first iteration. If None, the first
        dataset in `data` is used as reference. Default is None.
    n_iter : int, optional
        Number of iterations. Default is 10.
    tol : float, optional
        Tolerance for stopping criteria. Default is 1e-5.
    return_reference : bool, optional
        Whether to return the reference dataset built in the last iteration.
        Default is False.
    verbose : bool, optional
        Verbosity. Default is False.
    Returns
    -------
    aligned : list of ndarray, shape = (n_samples, n_feat)
        Aligned datsets.
    mean_dataset : ndarray, shape = (n_samples, n_feat)
        Reference dataset built in the last iteration. Only if
        ``return_reference == True``.
    """

    if n_iter <= 0:
        raise ValueError('A positive number of iterations is required.')

    if reference is None:
        # Use the first item to build the initial reference
        aligned = [data[0]] + [procrustes(d, data[0]) for d in data[1:]]
        reference = np.mean(aligned, axis=0)
    else:
        aligned = [None] * len(data)
        reference = reference.copy()

    dist = np.inf
    for i in range(n_iter):
        # Align to reference
        aligned = [procrustes(d, reference) for d in data]

        # Compute new mean
        new_reference = np.mean(aligned, axis=0)

        # Compute distance
        new_dist = np.square(reference - new_reference).sum()

        # Update reference
        reference = new_reference

        if verbose:
            print('Iteration {0:>3}: {1:.6f}'.format(i, new_dist))

        if dist != np.inf and np.abs(new_dist - dist) < tol:
            break

        dist = new_dist

    return (aligned, reference) if return_reference else aligned
  
def get_icc(measures):
    return icc.icc(measures,model='twoway',type='consistency',unit='single')
           
    
    