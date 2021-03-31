#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 12:34:21 2021

@author: bwinsto2
"""

import scipy.sparse
import scipy.sparse.linalg
import time
import numpy as np
import cupy

adj_mat = scipy.sparse.load_npz('/data/hcp_test_retest_pp/derivatives/chap/sub-103818/ses-test/struc_conn_mat.npz')
print(scipy.sparse.issparse(adj_mat))
print(adj_mat)

def lapDecomp(Asparse,num):
    #Asparse- csr format sparse adjacency matrix
    #num - number of laplacian eigenmodes to compute
    #returns tuple of eigenvalues (vals) and eigenvectors (vecs) as arrays of size (dim(Asparse),1) and  size (dim(Asparse),num)
    start = time.time()
    LA=scipy.sparse.csgraph.laplacian(Asparse,normed=True)
    L=LA.astype('float64')
    vals,vecs= scipy.sparse.linalg.eigsh(L,num,which='SM') #takes forever
    end = time.time()
    print('decomposition time=', (end-start),'seconds, which is', (end - start)/60, 'minutes')
    return vals,vecs
    
print(lapDecomp(adj_mat, 100))

#https://pyculib.readthedocs.io/en/latest/cusparse.html

