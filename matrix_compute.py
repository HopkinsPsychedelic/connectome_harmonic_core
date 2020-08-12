#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 18:27:26 2020

@author: patricktaylor
"""
from scipy import sparse
import numpy as np 
import utility_functions as ut
from sklearn.metrics import pairwise_distances_chunked
import time
from scipy.sparse import csgraph

def construct_surface_matrix(SC,SI):
    # SC- surface coordinates (used for determining size of surface matrix only, can be empty array of length=len(SC))
    # SI - array of vertex connections. Each row of SI contains indices of vertices in SC that form a triangle in the mesh
    M=sparse.lil_matrix((len(SC),len(SC))) #initialize
    #set each M_ij, where i,j are vertices connected in the mesh, equal to 1.
    M[SI[:,0],SI[:,1]]=1
    M[SI[:,1],SI[:,0]]=1
    M[SI[:,1],SI[:,2]]=1
    M[SI[:,2],SI[:,1]]=1
    M[SI[:,2],SI[:,0]]=1
    M[SI[:,0],SI[:,2]]=1
    return M.tocsr()


def construct_inter_hemi_matrix(SC,tol=4):
    '''
    creates interhemisphere connection matrix IHC for a given surface mesh with coordinates SC by connecting all vertices on the midline fissure 
    to their nearest neighbor on the opposite hemisphere. IHC has dimension (len(SC),len(SC)).
    '''
    half=int(len(SC)/2)
    li,ld=ut.neighbors(SC[:half],SC[half:],1)
    ri,rd=ut.neighbors(SC[half:],SC[:half],1)
    IHC=sparse.lil_matrix((half*2,half*2))
    #R=sparse.lil_matrix((half,half))
    for i in range (half):
        if ld[i]<tol:
            IHC[i+half,li[i]]=1
            IHC[li[i],i+half]=1
    for i in range (half):
        if rd[i]<tol:
            IHC[ri[i]+half,i]=1
            IHC[i,ri[i]+half]=1
    IHC=IHC.tocsr()
    return IHC



def construct_microstructure_matrix_inverse(msvecs, mc, SM, IHC,k=20):

    M = sparse.lil_matrix((len(mc), len(mc)))  # initialize empty result matrix
    r, c, v = sparse.find(SM + IHC)  # extract and flatten nonzero indices of SM
    dif = msvecs[r] - msvecs[c]
    dist = np.linalg.norm(dif, axis=1)
    for i in range(len(r)):

        M[r[i], c[i]] = (1 / (1 + k*dist[i]))

    x = np.arange(len(mc))
    M[x, x] = 0
    return M.tocsr()

def construct_structural_connectivity_matrix(SC,EC,tol=3,NNnum=45):
    ind,dist=ut.neighbors(SC,EC,1)
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

def construct_functional_connectivity_matrix_chunked(timeseries,threshold=0.01):
    starttime=time.time()
    #dense_covar=np.nan_to_num(np.corrcoef(timeseries))
    dense_covar=(np.corrcoef(timeseries))
    #flat_dense_covar=np.reshape(dense_covar,(len(dense_covar)*len(dense_covar)))
    #flat_dense_covar=np.nan_to_num(flat_dense_covar)
    #flat_zscore=scipy.stats.zscore(flat_dense_covar_nonan)
    #dense_covar_mat=np.reshape(flat_dense_covar,np.shape(dense_covar)) 
    #x=np.arange(len(dense_covar))
    #dense_covar[x,x]=0
    print('covariance done')
    sparse_threshold_covar=sparse.lil_matrix(np.shape(dense_covar))
    keep_num=int(np.round(len(dense_covar)*threshold))
    for i in range (len(dense_covar)):
        #print(i)
        row=np.nan_to_num(dense_covar[i,:])
        keep_inds=np.argpartition(row,-keep_num)[-keep_num:]
        sparse_threshold_covar[i,keep_inds]=row[keep_inds]
    
    sparse_threshold_covar=sparse_threshold_covar.tocsr()
    gen = pairwise_distances_chunked(sparse_threshold_covar,metric='cosine',n_jobs=-1)
    cosine_similarity=sparse.lil_matrix((len(timeseries),len(timeseries)))
    start=0
    for item in gen:
        end=start+len(item)
        g=1-item
        print(np.shape(g))
        cosine_similarity[start:end,:]=g
        start=end
    #cosine_similarity=1-cosine_distance
    x=np.arange(len(dense_covar))
    cosine_similarity[x,x]=0
    #cosine_similarity_sparse=sparse.csr_matrix(cosine_similarity)
    #symmetric_sparse_threshold_covar=(sparse_threshold_covar+sparse_threshold_covar.T)/2
    endtime=time.time()
    print(endtime-starttime)
    return cosine_similarity.tocsr()

def construct_functional_connectivity_matrix_chunked_transpose(timeseries,threshold=0.01):
    starttime=time.time()
    #dense_covar=np.nan_to_num(np.corrcoef(timeseries))
    dense_covar=(np.corrcoef(timeseries))
    #flat_dense_covar=np.reshape(dense_covar,(len(dense_covar)*len(dense_covar)))
    #flat_dense_covar=np.nan_to_num(flat_dense_covar)
    #flat_zscore=scipy.stats.zscore(flat_dense_covar_nonan)
    #dense_covar_mat=np.reshape(flat_dense_covar,np.shape(dense_covar)) 
    #x=np.arange(len(dense_covar))
    #dense_covar[x,x]=0
    print('covariance done')
    sparse_threshold_covar=sparse.lil_matrix(np.shape(dense_covar))
    keep_num=int(np.round(len(dense_covar)*threshold))
    for i in range (len(dense_covar)):
        #print(i)
        row=np.nan_to_num(dense_covar[i,:])
        keep_inds=np.argpartition(row,-keep_num)[-keep_num:]
        sparse_threshold_covar[i,keep_inds]=row[keep_inds]
    
    sparse_threshold_covar=sparse_threshold_covar.tocsr()
    
    x=np.arange(np.shape(sparse_threshold_covar)[0])
    sparse_threshold_covar[x,x]=0
    
    #cosine_similarity_sparse=sparse.csr_matrix(cosine_similarity)
    #symmetric_sparse_threshold_covar=(sparse_threshold_covar+sparse_threshold_covar.T)/2
    endtime=time.time()
    print(endtime-starttime)
    return (sparse_threshold_covar+sparse_threshold_covar.T)/2

def construct_single_scalar_matrix_inverse(scalar, SM, ihc):
    M = sparse.lil_matrix((len(scalar), len(scalar)))
    r, c, v = sparse.find(SM + ihc)
    dif = scalar[r] - scalar[c]
    dot = np.square(dif)
    dist = np.sqrt(dot)
    for i in range(len(r)):
        M[r[i], c[i]] = (1 / (1 + dist[i]))
    M.data = np.nan_to_num(M.data)
    x = np.arange(len(scalar))
    M[x, x] = 0
    return M.tocsr()

def construct_microstructure_matrix_gaussian(msvecs,SM,IHC,kappa=.5): 
    '''
    msvecs- microstructure vertex-wise feature vectors
    mc- surface coordinates 
    SM- surface matrix used to define local neighborhood of each vertex
    IHC- interhemicon matrix, same use as SM
    kappa- used in kernelized weight calculation: w_ij=C*exp(kappa*(<v_i,v_j>)^2), 
        where w_ij is the connection weight and <v_i,v_j> is the cosine norm between the feature vectors of vertex i and vertex j
    
    '''
    M=sparse.lil_matrix((len(msvecs),len(msvecs))) #initialize empty result matrix
    r,c,v=sparse.find(SM) #extract and flatten nonzero indices of SM
    a=np.einsum('ij, ij->i', msvecs[r], msvecs[r]) # compute magnitude of feature vectors
    b=np.einsum('ij, ij->i', msvecs[c], msvecs[c]) # ^^^^      ^^^^     ^^
    d=np.einsum('ij, ij->i',msvecs[r] ,msvecs[c])  # compute unormalized dot product between all vertices connected by SM
    M[r,c]=np.nan_to_num((1/np.exp(kappa))*np.exp(-kappa*np.square(np.true_divide(d,np.multiply(np.sqrt(a),np.sqrt(b)))))) #populate result matrix M using gaussian kernel of cosine norm
    R,C,V=sparse.find(IHC) #repeat above process with IHC connection matrix
    A=np.einsum('ij, ij->i', msvecs[R], msvecs[R])
    B=np.einsum('ij, ij->i', msvecs[C], msvecs[C])
    D=np.einsum('ij, ij->i',msvecs[R] ,msvecs[C])
    M[R,C]=np.nan_to_num((1/np.exp(kappa))*np.exp(-kappa*np.square(np.true_divide(D,np.multiply(np.sqrt(A),np.sqrt(B))))))
    
    M.data=np.nan_to_num(M.data)
    M=np.nan_to_num(M)
    rr,cc,vv=sparse.find(M)
    print('mean val=',np.mean(vv),'median val=',np.median(vv))
    print('max val=',np.max(vv),'min val=',np.min(vv)) 
    M=M.tocsr()
    
    return M

def diffusion_matrix(A,t=0):
    Lap,D_sqrt_vec=csgraph.laplacian(A,normed=True,return_diag=True)
    D_inv_sqrt_vec=1/D_sqrt_vec
    D_inv_sqrt=sparse.lil_matrix(np.shape(Lap))
    x=np.arange(np.shape(Lap)[0])
    D_inv_sqrt[x,x]=D_inv_sqrt_vec
    L_sqrt=(D_inv_sqrt.dot(A)).dot(D_inv_sqrt)
    Lap,D_alpha=csgraph.laplacian(L_sqrt,return_diag=True)
    D_alpha_mat=sparse.lil_matrix(np.shape(Lap))
    D_alpha_mat[x,x]=1/D_alpha
    P=D_alpha_mat.dot(L_sqrt)
    print('diffusion matrix computed')
    if t==0:
        return P
    else:
        return P**t

def cosine_similarity_matrix(M):
    gen = pairwise_distances_chunked(M,metric='cosine',n_jobs=-1)
    cosine_similarity=sparse.lil_matrix(np.shape(M))
    start=0
    for item in gen:
        end=start+len(item)
        g=1-item
        #g=np.where(g<0.1,0,g)
        print(np.shape(g))
        cosine_similarity[start:end,:]=g
        start=end
    #cosine_similarity=1-cosine_distance
    x=np.arange(np.shape(M)[0])
    cosine_similarity[x,x]=0
    return cosine_similarity

