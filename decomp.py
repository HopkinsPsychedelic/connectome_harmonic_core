#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 19:36:03 2020

@author: patricktaylor
"""
import scipy.sparse.linalg
from scipy.sparse import csgraph
import time
import numpy as np

from sklearn.decomposition import PCA

def get_group_pca_comp(evlist,num):
    tempmat = np.zeros((np.shape(evlist[0])[1]-1,num*len(evlist)))
    for i in range (len(evlist)):
        tempmat[:,i*num:(i+1)*num]=evlist[i][:,1:num] #from 0-10, 10-20, 20-30   
    pca = PCA(n_components=num)
    pca.fit(tempmat.T)
    components = pca.components_
    components = components.T
    return components 

def get_group_pca_comp_b(evlist,num):
    tempmat = np.zeros((np.shape(evlist[0])[0],num*len(evlist)))
    for i in range (len(evlist)):
        tempmat[:,i*num:(i+1)*num]=evlist[i][:,0:num] #from 0-10, 10-20, 20-30   
    pca = PCA(n_components=num)
    pca.fit(tempmat.T)
    components = pca.components_
    components = components.T
    return components 

def lapDecomp(Asparse,num):
    #Asparse- csr format sparse adjacency matrix
    #num - number of laplacian eigenmodes to compute
    #returns tuple of eigenvalues (vals) and eigenvectors (vecs) as arrays of size (dim(Asparse),1) and  size (dim(Asparse),num)
    start = time.time()
    LA=csgraph.laplacian(Asparse,normed=True)
    L=LA.astype('float64')
    vals,vecs=scipy.sparse.linalg.eigsh(L,num,which='SM')
    end = time.time()
    print('decomposition time=', (end-start),'seconds, which is', (end - start)/60, 'minutes')
    return vals,vecs
