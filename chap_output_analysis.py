#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 09:46:54 2020

@author: bwinston
"""
import numpy as np
import os
from glob import glob
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
import statistics as stats     
         
'''
matching algorithm
for every ses 1 harmonic, correlation with all others
standard deviation and error, all harms with .8, .7, sequentially
each subj. come up with diagonal of best pairs 
components across subjects (get set of reliable harmonics, look across subjects)
dont steal from winners -- tiebreaker can be filled by another tie or loser 

average two sessions components? maybe
do both average of two sessions and just first session and compare to group average connectome harmonics
continuous variable on non binary matrix - how to do laplacian for average connectome 
get rid of noise components--those not reliable within a subject or across subjects. interesting will be something
that's reliable within a subject but not across subjects. (that's where individual differences lie)
so how do we find reliable harmonics across subjects?
then we'll find test retest reliability of that, which is a more accurate depiction of test retest rel
fxnl - 
What would you be plugging in to PCA?
PCA would be on set of all harmonics 


'''

test_retest_rel('/Users/bwinston/Documents/connectome_harmonics/chap_output/chap', 200) 
test_retest_rel('/Users/bwinston/Downloads/chap_out_test', 200)    
 
def get_key(my_dict, val):
    for key, value in my_dict.items():
         if val == value:
             return key
         
def test_retest_rel(chap_dir, n_evecs):
    global cd 
    cd = {} #for chap_data
    all_bcorrs, bcorr_plot = [], []
    subject_dirs = glob(os.path.join(chap_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
    for sub in subs:
        cd[sub] = {}
        cd[sub]['bcorrs'] = []
        for ses in ['test','retest']:
           cd[sub][ses] = {}
           cd[sub][ses]['vecs'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/vecs.npy') 
        hp, hp['holy'] = {}, {}
        hp['corr_orig'] = np.empty((n_evecs,n_evecs))
        for evec_test in range(0,n_evecs): 
            for evec_retest in range(0,n_evecs): #n_evecs x n_evecs correlation matrix
                hp['corr_orig'][evec_retest,evec_test] = abs(pearsonr(cd[sub]['test']['vecs'][:,evec_test], cd[sub]['retest']['vecs'][:,evec_retest])[0])
        hp['corr_all'] = hp['corr_orig'].swapaxes(0,1)
        hp['corr_all'] = {index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(hp['corr_all'])}
        find_bcorrs(sub, hp, 0, hp['corr_all'], n_evecs)
        for ev in range(n_evecs):
            cd[sub]['bcorrs'].append(cd[sub]['pairs'][ev]['bcorr'])
    for sub in subs:        
        all_bcorrs.append(cd[sub]['bcorrs'])
    cd['bcorr_avg'] = np.average(np.array(all_bcorrs), axis=0)
    for avg in range(1,n_evecs):
        bcorr_plot.append(stats.mean(cd['bcorr_avg'][:avg])) 
    plt.plot(bcorr_plot)

def find_bcorrs(sub, hp, run, corr_mat, n_evecs): 
    while len(hp['corr_all']) > 0:
        hp[run] = {}
        hp[run]['maxes'] = {}
        for ev in hp['corr_all']: #for each test session evec, which retest session evec is the best match?
            hp[run]['maxes'][ev] = {}
            hp[run]['maxes'][ev]['bcorr'] = max(hp['corr_all'][ev].values())
            hp[run]['maxes'][ev]['ret_ind'] = get_key(hp['corr_all'][ev], max(hp['corr_all'][ev].values()))
        hp[run]['ret_used'], hp[run]['win_ind'], hp[run]['good_boys'] = [], [], []
        for ev in hp['corr_all']:
            hp[run]['ret_used'].append(hp[run]['maxes'][ev]['ret_ind'])
        ret_dups = set([x for x in hp[run]['ret_used'] if hp[run]['ret_used'].count(x) > 1])
        if ret_dups == 0:
            break
        hp[run]['ret_used'] = set(hp[run]['ret_used'])
        for ret_dup in ret_dups:
            competition = []
            for ev in hp[run]['maxes']:
                if ret_dup == hp[run]['maxes'][ev]['ret_ind']:
                    competition.append(hp[run]['maxes'][ev]['bcorr'])
            competition = max(competition)
            for ev in hp[run]['maxes']:
                if hp[run]['maxes'][ev]['bcorr'] == competition:
                    hp[run]['win_ind'].append(ev)
        for ev in hp['corr_all']:
            if hp[run]['maxes'][ev]['ret_ind'] not in ret_dups:
               hp[run]['good_boys'].append(ev) 
        for ev in list(hp['corr_all']):
            if ev in (set(hp[run]['good_boys']) | set(hp[run]['win_ind'])):
               hp['holy'][ev] = hp[run]['maxes'][ev]
               del hp['corr_all'][ev]
        for ev in hp['corr_all']:
            for r_ev in hp[run]['ret_used']:
                del hp['corr_all'][ev][r_ev]
        run = run + 1
        find_bcorrs(sub, hp, run, hp['corr_all'], len(hp['corr_all']))
    cd[sub]['pairs'] = hp['holy']
    
'''
average shit
'''

def avg_harms(chap_dir):
    subject_dirs = glob(os.path.join(chap_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]
    global avg
    avg = {}
    avg['connectomes_test'], avg['connectomes_retest'], avg['connectomes_all']= 0,0,0
    for sub in subs:
       avg['connectomes_test'] = avg['connectomes_test'] + sparse.load_npz(f'{chap_dir}/sub-{sub}/ses-test/connectome.npz') #sum all test connectomes
       avg['connectomes_retest'] = avg['connectomes_retest'] + sparse.load_npz(f'{chap_dir}/sub-{sub}/ses-retest/connectome.npz')
       avg['connectomes_all'] = avg['connectomes_all'] + sparse.load_npz(f'{chap_dir}/sub-{sub}/ses-test/connectome.npz') + sparse.load_npz(f'{chap_dir}/sub-{sub}/ses-retest/connectome.npz')
    avg['connectomes_test'] = avg['connectomes_test'] / len(subs) #take average
    avg['connectomes_retest'] = avg['connectomes_retest'] / len(subs)
    avg['connectomes_all'] = avg['connectomes_all'] / (len(subs)*2)
    inout.if_not_exist_make(f'{chap_dir}/sub-test_avg')
    inout.if_not_exist_make(f'{chap_dir}/sub-retest_avg')
    inout.if_not_exist_make(f'{chap_dir}/sub-total_avg')
    avg_test_vals,avg_test_vecs = dcp.lapDecomp(avg['connectomes_test'], 200) 
    avg_retest_vals,avg_retest_vecs = dcp.lapDecomp(avg['connectomes_retest'], 200) 
    avg_total_vals,avg_total_vecs = dcp.lapDecomp(avg['connectomes_all'], 200) 
    np.save(f'{chap_dir}/sub-test_avg/vecs', avg_test_vecs)
    np.save(f'{chap_dir}/sub-retest_avg/vecs', avg_retest_vecs)
    np.save(f'{chap_dir}/sub-total_avg/vecs', avg_total_vecs)

'''

 
'''
def ind_vs_avg(chap_dir, n_evecs):
    global iva
    iva = {} 
    subject_dirs = glob(os.path.join(chap_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs]
    for sub in ['test_avg', 'retest_avg', 'total_avg']:
        subs.remove(sub)
    for sub in subs:
        iva[sub] = {}    
        for ses in ['test','retest']:
           iva[sub][ses],iva[f'{ses}_avg'] = {}, {}
           iva[sub][ses]['bcorrs'] = []
           iva[sub][ses]['vecs'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/vecs.npy')
           iva[f'{ses}_avg']['vecs'] = np.load(f'{chap_dir}/sub-{ses}_avg/vecs.npy')
        global hp
        hp = {}
        for ses in ['test','retest']:
            hp[ses], hp[ses]['holy'] = {}, {}  
            hp[ses]['corr_orig'] = np.empty((n_evecs,n_evecs)) #init comparison matrix 
            for evec_avg in range(0,n_evecs): 
                for evec_test in range(0,n_evecs): #n_evecs x n_evecs correlation matrix
                    hp[ses]['corr_orig'][evec_test,evec_avg] = abs(pearsonr(iva['test_avg']['vecs'][:,evec_avg], iva[sub][ses]['vecs'][:,evec_test])[0])
            hp[ses]['corr_all'] = hp[ses]['corr_orig'].swapaxes(0,1)
            hp[ses]['corr_all'] = {index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(hp[ses]['corr_all'])}
            find_bcorrs_iva(sub, ses, hp, 0, hp[ses]['corr_all'], n_evecs)
    for sub in subs:
        for ses in ['test','retest']:
            for ev in range(n_evecs):
                iva[sub][ses]['bcorrs'].append(iva[sub][ses]['pairs'][ev]['bcorr'])
    iva['all_test_bcorrs'],iva['all_retest_bcorrs'] = [],[]
    for sub in subs:
        for ses in ['test','retest']:
            iva[f'all_{ses}_bcorrs'].append(iva[sub][ses]['bcorrs'])
    iva['bcorr_test_avg'] = np.average(np.array(iva['all_test_bcorrs']), axis = 0)
    iva['bcorr_retest_avg'] = np.average(np.array(iva['all_retest_bcorrs']), axis = 0)
        
        

def find_bcorrs_iva(sub, ses, hp, run, corr_mat, n_evecs): 
    while len(hp[ses]['corr_all']) > 0:
        hp[ses][run] = {}
        hp[ses][run]['maxes'] = {}
        for ev in hp[ses]['corr_all']: #for each test session evec, which retest session evec is the best match?
            hp[ses][run]['maxes'][ev] = {}
            hp[ses][run]['maxes'][ev]['bcorr'] = max(hp[ses]['corr_all'][ev].values())
            hp[ses][run]['maxes'][ev][f'{ses}_ind'] = get_key(hp[ses]['corr_all'][ev], max(hp[ses]['corr_all'][ev].values()))
        hp[ses][run][f'{ses}_used'], hp[ses][run]['win_ind'], hp[ses][run]['good_boys'] = [], [], []
        for ev in hp[ses]['corr_all']:
            hp[ses][run][f'{ses}_used'].append(hp[ses][run]['maxes'][ev][f'{ses}_ind'])
        ses_dups = set([x for x in hp[ses][run][f'{ses}_used'] if hp[ses][run][f'{ses}_used'].count(x) > 1])
        if ses_dups == 0:
            break
        hp[ses][run][f'{ses}_used'] = set(hp[ses][run][f'{ses}_used'])
        for ses_dup in ses_dups:
            competition = []
            for ev in hp[ses][run]['maxes']:
                if ses_dup == hp[ses][run]['maxes'][ev][f'{ses}_ind']:
                    competition.append(hp[ses][run]['maxes'][ev]['bcorr'])
            competition = max(competition)
            for ev in hp[ses][run]['maxes']:
                if hp[ses][run]['maxes'][ev]['bcorr'] == competition:
                    hp[ses][run]['win_ind'].append(ev)
        for ev in hp[ses]['corr_all']:
            if hp[ses][run]['maxes'][ev][f'{ses}_ind'] not in ses_dups:
               hp[ses][run]['good_boys'].append(ev) 
        for ev in list(hp[ses]['corr_all']):
            if ev in (set(hp[ses][run]['good_boys']) | set(hp[ses][run]['win_ind'])):
               hp[ses]['holy'][ev] = hp[ses][run]['maxes'][ev]
               del hp[ses]['corr_all'][ev]
        for ev in hp[ses]['corr_all']:
            for s_ev in hp[ses][run][f'{ses}_used']:
                del hp[ses]['corr_all'][ev][s_ev]
        run = run + 1
        find_bcorrs_iva(sub, ses, hp, run, hp[ses]['corr_all'], len(hp[ses]['corr_all']))
    iva[sub][ses]['pairs'] = hp[ses]['holy']



po = iva['bcorr_test_avg']
op = iva['bcorr_retest_avg']

po = list(po)
po.sort(reverse=True)

avg_harms('/Users/bwinston/Documents/connectome_harmonics/chap_output/chap')
ind_vs_avg('/Users/bwinston/Documents/connectome_harmonics/chap_output/chap', 100)


#ask patrick if visualization is possible
#can make an average struc conn mat, but having an average sc and si doesn't really make sense

   inout.save_eigenvector(f'{args.output_dir}/chap/sub-{sub}/{ses}/vis/sub-{sub}_{ses}_harmonics.vtk',sc,si,vecs) #harmonics.vtk
        





hi = np.load('/Users/bwinston/Documents/connectome_harmonics/chap_output/avg_stuff/avg_test_vecs.npy')
np.shape(hi[:,0])
yo = np.load('/Users/bwinston/Documents/connectome_harmonics/chap_output/chap/sub-105923/ses-test/vecs.npy')
np.shape(yo[:,0])
type(hi[:,0])







     
