#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 13:00:45 2020

@author: bwinston
"""
import numpy as np
import matrix_methods as mm
import input_output as inout
import decomp as dcp
from scipy import sparse
import numpy as np
import utility_functions as ut
from sklearn.metrics import pairwise_distances_chunked
import time
from scipy.sparse import csgraph
import matplotlib.pylab as plt
from scipy.stats.stats import pearsonr
import statistics

vecs = np.load('/Users/bwinston/Documents/connectome_harmonics/vecs.npy')
vals = np.load('/Users/bwinston/Documents/connectome_harmonics/vals.npy')

sc,si = inout.read_gifti_surface_both_hem('/Users/bwinston/Documents/fMRI/BIDS/HCP_Preproc/ses-retest/105923/T1w/fsaverage_LR32k/105923.L.white.32k_fs_LR.surf.gii', '/Users/bwinston/Documents/fMRI/BIDS/HCP_Preproc/ses-retest/105923/T1w/fsaverage_LR32k/105923.R.white.32k_fs_LR.surf.gii', hcp=True)

ec = inout.read_streamline_endpoints('/Users/bwinston/Downloads/10000000_endpoints.vtk')

def brian_sc_mat(SC,EC,tol=3,NNnum=45):
    '''
    SC- array of cortical surface coordinates of size (N_vertices, 3 ) where SC[i]=x_i,y_i,z_i
    EC- array of streamline endpoint coordinates of size (2*N_streamlines, 3 ) where EC[i]=[x_i,y_i,z_i]. 
    tol- search radius of nearest neighbor search for matching endpoints to surface vertices
    NNnum- number of nearest neighboring surface vertices to assign to each endpoint
    '''
    ind,dist=ut.neighbors(SC,EC,1)
    print(f'ind len is {ind.shape} and dist len is {dist.shape}')
    print(f'ind first 10 is {ind[:10]}; dist first ten is {dist[:10]}')
    bad=[]
    c=np.arange(len(dist)) 
    even=c[::2] 
    odd=c[1::2]
    for i in range (int(len(dist)/2)): 
        if (dist[even[i]]>=tol or dist[odd[i]]>=tol): 
            bad.append(even[i])
            bad.append(odd[i])
    newEC=np.delete(EC,bad,axis=0)
    s2eInd, s2eDist=ut.neighbors(SC,newEC,1)
    Rind,Rdist=ut.neighbors(newEC,SC,NNnum)
    OtherEndInd=np.zeros(np.shape(Rind))
    for i in range(len(Rind)):
        for j in range (NNnum):
            if Rind[i][j]%2==0 :
                OtherEndInd[i][j]=int(s2eInd[Rind[i][j]+1])
            else:
                OtherEndInd[i][j]=int(s2eInd[Rind[i][j]-1])

    M=sparse.lil_matrix((len(SC),len(SC)))
    x=np.arange(len(SC))
    for i in range (NNnum):
        AccSurfInd=np.column_stack((x,OtherEndInd[:,i]))
        U,C=np.unique(AccSurfInd,axis=0,return_counts=True)
        M[U[:,0],U[:,1]]+=C
        M[U[:,1],U[:,0]]+=C
    print(M.nnz)

    x=np.arange(len(SC))
    M[x,x]=0
    print(M.nnz)
    M=M.tocsr()
    return M.tocsr()

struc = brian_sc_mat(sc, ec)
plt.spy(struc, markersize = .005)
surf_mat=mm.construct_surface_matrix(sc,si)
plt.spy(surf_mat, markersize = .005)

connectome = struc + surf_mat

hi = connectome.todense()

plt.spy(connectome, markersize = .005)


hi = struc.tolil()
print(hi[4])

x = np.arange(1000)
even = x[::2]
odd = x[1::2]

#pearson correlations btwn vectors
vecs_test = np.load('/Users/bwinston/Downloads/sub-122317/ses-test/vecs.npy')
vecs_retest = np.load('/Users/bwinston/Downloads/sub-122317/ses-retest/vecs.npy')
vecs_test[:,0] #gets column or evec
corr_200109 = np.empty((200,200))
for evec_test in range(0,200): 
    for evec_retest in range(0,200): #200x200 correlation of each evec w/ each other evec
        corr_200109[evec_retest,evec_test] = abs(pearsonr(vecs_test[:,evec_test], vecs_retest[:,evec_retest])[0])
plt.imshow(corr_200109, cmap= 'plasma')
plt.colorbar()
plt.show()
self_harms = np.empty((1,200))
for corr in range(0,200): #correlation btwn test_evec_i and retest_evec_i from 1-200
    self_harms[:,corr] = corr_200109[corr][corr]

mean_of_harm_corr = []
for x in range(1,200):
    mean_of_harm_corr.append(statistics.mean(self_harms[0][:x]))
plt.plot(mean_of_harm_corr, 'b')

test_maxes = {}
test_maxes['corr'],  test_maxes['ind'] = [], [] #for each test session evec, which retest session evec is the best match?
for corr in range(len(corr_200109[0])):
    test_maxes['corr'].append(max(corr_200109[:,corr])) #maximum correlation of that evec to another evec
    test_maxes['ind'].append(np.where(corr_200109[:,corr]==test_maxes['corr'][corr])[0][0]) #which evec from other session does the max correspond to

retest_maxes = {}
retest_maxes['corr'],  retest_maxes['ind'] = [], [] #for each retest session evec, which test session evec is the best match?
for corr in range(len(corr_200109[0])):
    retest_maxes['corr'].append(max(corr_200109[corr])) #maximum correlation of that evec to another evec
    retest_maxes['ind'].append(np.where(corr_200109[corr]==retest_maxes['corr'][corr])[0][0]) #which evec from other session does the max correspond to

mean_of_good_fits = []
for x in range(1,200):
    mean_of_good_fits.append(statistics.mean(retest_maxes['corr'][:x]))
plt.plot(mean_of_good_fits, 'g')


        
        
        
    

