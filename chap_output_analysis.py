#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 09:46:54 2020

@author: bwinston
"""
import numpy as np
import os
from glob import glob
#import matrix_methods as mm
#import input_output as inout
#import decomp as dcp
from scipy import sparse
import numpy as np
#import utility_functions as ut
from sklearn.metrics import pairwise_distances_chunked
import time
from scipy.sparse import csgraph
import matplotlib.pylab as plt
from scipy.stats.stats import pearsonr
import statistics as stats


def test_retest_rel(chap_dir, num_evecs):
    global cd 
    cd = {} #for chap_data
    subject_dirs = glob(os.path.join(chap_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
    cd['corr_all'], bcorrs_test_avg, bcorrs_retest_avg, cd['all_bcorrs_test'], cd['all_bcorrs_retest'] = [], [], [], [], []
    for sub in subs:
        cd[sub] = {}
        for ses in ['test','retest']:
           cd[sub][ses] = {}
           cd[sub][ses]['vecs'] = np.load(f'{chap_dir}/sub-{sub}/ses-{ses}/vecs.npy')
        cd[sub]['corr'] = np.empty((num_evecs, num_evecs))
        for evec_test in range(0,num_evecs):
            for evec_retest in range(0,num_evecs):
                cd[sub]['corr'][evec_retest,evec_test] = abs(pearsonr(cd[sub]['test']['vecs'][:,evec_test], cd[sub]['retest']['vecs'][:,evec_retest])[0])
        plt.imshow(cd[sub]['corr'], cmap = 'cividis')
        plt.colorbar()
        plt.show()
        cd['corr_all'].append(cd[sub]['corr'])
        cd[sub]['self_harms'], cd[sub]['mean_of_harm_conn']  = [], []
        for i in range(0,num_evecs):
            cd[sub]['self_harms'].append(cd[sub]['corr'][i][i])
            if i>0:
                cd[sub]['mean_of_harm_conn'].append(stats.mean(cd[sub]['self_harms']))
        for ses in ['test','retest']:
            cd[sub][ses]['maxes'] = {}
            cd[sub][ses]['maxes']['bcorr'], cd[sub][ses]['maxes']['ind'] = [], []
        for ses in ['test', 'retest']:      
            for corr in range(0, num_evecs):
                if ses == 'test':
                    cd[sub][ses]['maxes']['bcorr'].append(max(cd[sub]['corr'][:,corr]))
                else:
                    cd[sub][ses]['maxes']['bcorr'].append(max(cd[sub]['corr'][corr]))
                cd[sub][ses]['maxes']['ind'].append(np.where(cd[sub]['corr'][:,corr]==cd[sub][ses]['maxes']['bcorr'][corr])[0])
                cd[sub][ses]['maxes']['mogf'] = []
                for nharms in range(1,num_evecs):
                    cd[sub][ses]['maxes']['mogf'].append(stats.mean(cd[sub][ses]['maxes']['bcorr'][:nharms]))
            cd[f'all_bcorrs_{ses}'].append(cd[sub][ses]['maxes']['bcorr'])
    cd['corr_total_avg'] = np.mean(np.array(cd['corr_all']), axis=0)
    plt.imshow(cd['corr_total_avg'], cmap = 'plasma')
    plt.colorbar()
    plt.show()
    cd['all_bcorrs_test'] = np.array(cd['all_bcorrs_test'])
    cd['all_bcorrs_retest'] = np.array(cd['all_bcorrs_retest'])
    cd['bcorr_test_avg'] = np.mean(np.array(cd['all_bcorrs_test']), axis=0)
    cd['bcorr_retest_avg'] = np.mean(np.array(cd['all_bcorrs_retest']), axis=0)
    for harm in range(1,num_evecs):
        bcorrs_test_avg.append(stats.mean(cd['bcorr_test_avg'][:harm]))   
        bcorrs_retest_avg.append(stats.mean(cd['bcorr_retest_avg'][:harm]))
    plt.plot(bcorrs_test_avg)
    plt.plot(bcorrs_retest_avg)
    
test_retest_rel('/Users/bwinston/Downloads/chap_out_test', 200)       
         
'''
matching algorithm
for every ses 1 harmonic, correlation with all others
standard deviation and error, all harms with .8, .7, sequentially
each subj. come up with diagonal of best pairs 
components across subjects (get set of reliable harmonics, look across subjects)
dont steal from winners -- tiebreaker can be filled by another tie or loser 

'''



sub = str(200109)
n_evecs = 200
vecs_test = np.load('/Users/bwinston/Downloads/chap_out_test/sub-200109/ses-test/vecs.npy')
vecs_retest = np.load('/Users/bwinston/Downloads/chap_out_test/sub-200109/ses-retest/vecs.npy')
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
    mean_of_harm_corr.append(stats.mean(self_harms[0][:x]))
plt.plot(mean_of_harm_corr, 'b')

test_maxes = {}
test_maxes['corr'],  test_maxes['ind'] = [], [] #for each test session evec, which retest session evec is the best match?
for corr in range(len(corr_200109[0])):
    test_maxes['corr'].append(max(corr_200109[:,corr])) #maximum correlation of that evec to another evec
    test_maxes['ind'].append(np.where(corr_200109[:,corr]==test_maxes['corr'][corr])[0][0]) #which evec from other session does the max correspond to

def get_key(my_dict, val):
    for key, value in my_dict.items():
         if val == value:
             return key

'''
important shit
'''

test_retest_rel('/Users/bwinston/Documents/connectome_harmonics/chap_output/chap', 1000)

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
        
        
for ev in range(200):  
    print(cd['122317']['pairs'][ev]['bcorr']) 

cd['122317']['pairs'][199]     
    

 

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

     
    
        



         


mean_of_good_fits = []
for x in range(1,200):
    mean_of_good_fits.append(stats.mean(retest_maxes['corr'][:x]))
plt.plot(mean_of_good_fits, 'g')
