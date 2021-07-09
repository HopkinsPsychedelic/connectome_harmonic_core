import numpy as np
from scipy import sparse 
import time 

def normalize_time_series(ts):
    ts=ts.T
    nts=(ts - np.mean(ts, axis=0)) / np.std(ts, axis=0)
    return nts.T

def truncate_top_k(x, k, inplace=False):
    m, n = x.shape
    # get (unsorted) indices of top-k values
    topk_indices = np.argpartition(x, -k, axis=1)[:, -k:]
    # get k-th value
    rows, _ = np.indices((m, k))
    kth_vals = x[rows, topk_indices].min(axis=1)
    # get boolean mask of values smaller than k-th
    is_smaller_than_kth = x < kth_vals[:, None]
    # replace mask by 0
    if not inplace:
        return np.where(is_smaller_than_kth, 0, x)
    x[is_smaller_than_kth] = 0
    return x

def construct_FC_matrix_row_thresh(ts,threshold=0.01,chunk_size=1000,normalize=True):
    
    start=time.time()
    
    if normalize:
        nts = normalize_time_series(ts).T
    else:
        nts =ts.T

    sparse_chunck_list = []
    
    keep_num=int(np.round(len(ts)*threshold))
    print(keep_num, 'elements per row retained')
    
    for i in range(0, int(nts.shape[1]), chunk_size):
        # portion of connectivity
        pcon = (np.matmul(nts[:, i:i + chunk_size].T, nts) / nts.shape[0])
        
        # sparsified connectivity portion
        spcon = sparse.csr_matrix(truncate_top_k(pcon,keep_num,inplace=True))
        
        
        
        sparse_chunck_list.append(spcon)
        
        print(f'rows {i}:{i+chunk_size} done')
        
    scon = sparse.vstack(sparse_chunck_list)
    
    #scon=scon.tolil()
    #scon[scon<0]=0
    scon=scon.tocsr()
    scon=(scon+scon.T)/2
    end = time.time()
    
    x=np.arange(len(ts))
    scon[x,x]=0
    print(f'{end-start} seconds taken')
    return scon.astype('float32')