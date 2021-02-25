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
import test_retest_fxns as t_rt
import utility_functions as ut
from sklearn.metrics import pairwise_distances_chunked, pairwise
import time
import matplotlib.pylab as plt
from scipy.stats.stats import pearsonr
from scipy.spatial import distance
import statistics as stats     
import compute_spectra as cs
         
'''
standard deviation and error, all harms with .8, .7, sequentially 
components across subjects (get set of reliable harmonics, look across subjects) 

average two sessions components? maybe
do both average of two sessions and just first session and compare to group average connectome harmonics
get rid of noise components--those not reliable within a subject or across subjects. interesting will be something
that's reliable within a subject but not across subjects. (that's where individual differences lie)
so how do we find reliable harmonics across subjects?
then we'll find test retest reliability of that, which is a more accurate depiction of test retest rel
PCA would be on set of all harmonics 


WITHIN SUBJECT TEST-RETEST RELIABILITY:
'''
t_rt.test_retest_rel('/Users/bwinston/Documents/connectome_harmonics/chap_output/chap', 200) 
t_rt.test_retest_rel('/Users/bwinston/Downloads/chap_out_test', 50) 
t_rt.test_retest_rel_r('/data2/Brian/connectome_harmonics/three_chap_subjs', 100)
t_rt.test_retest_rel_cos('/data2/Brian/connectome_harmonics/three_chap_subjs', 100)



    

hi = iva['105923']['test']['pairs']

po = iva['bcorr_test_avg']
op = iva['bcorr_retest_avg']

hcp_fem_ids = [103818, 105923, 111312, 114823, 115320, 125525, 130518, 135528, 137128, 143325, 144226, 158035, 169343, 172332, 175439, 177746, 187547, 192439, 194140, 195041, 200109, 200614, 204521, 250427, 287248, 562345,627549, 660951, 859671,861456, 877168]
hcp_male_ids = [122317, 139839, 146129, 149337, 149741, 151526, 185442, 341834, 433839, 599671, 601127, 783462, 917255]
hcp_all = hcp_fem_ids + hcp_male_ids
hcp_all = [str(i) for i in hcp_all]

avg_harms('/Users/bwinston/Documents/connectome_harmonics/chap_output/chap')
t_rt.ind_vs_avg('/Users/bwinston/Documents/connectome_harmonics/chap_output/chap', 100)

   
t_rt.ind_vs_pca('/Users/bwinston/Downloads/chap_out_test', 200, 100)   
t_rt.ind_vs_pca('/data/hcp_test_retest_pp/derivatives/chap', 200, 100)   
    
    
'''cosine similarity'''
 
vec_1 = np.load('/data2/Brian/connectome_harmonics/three_chap_subjs/sub-187547/ses-test/vecs.npy')[:,5]
vec_2 = np.load('/data2/Brian/connectome_harmonics/three_chap_subjs/sub-187547/ses-retest/vecs.npy')[:,5]
vec_1=vec_1.reshape(-1,1)
vec_2=vec_2.reshape(-1,1)
vec_1 = vec_1.swapaxes(0,1)
vec_2 = vec_2.swapaxes(0,1)
pairwise.cosine_similarity(vec_1, vec_2)    
def chap_cosine_sim(vec_1, vec_2):
    vec_1 = vec_1.reshape(-1,1)
    vec_1 = vec_1.swapaxes(0,1)
    vec_2 = vec_2.reshape(-1,1)
    vec_2 = vec_2.swapaxes(0,1)
    return abs(pairwise.cosine_similarity(vec_1, vec_2))        
chap_cosine_sim(vec_1, vec_2)   

'''Measures of structural test-retest reliability'''
#best correlations within subj. vs. across subj.
#significance of difference?
metric_1_forreal = t_rt.struc_metric_1('/data/hcp_test_retest_pp/derivatives/chap', 100)
metric_1 = np.array([metric_1_forreal[0], metric_1_forreal[1]])
np.save('/data/hcp_test_retest_pp/derivatives/chap_analysis/struc_metrics/struc_metric_1', metric_1)
#within subj. avg = .45; across subj. average = .25 (for three subjs, 100 evecs)
#both numbers get higher when you include fewer evecs

#where does within subj vs. across subj. separate?
sep_forreal = t_rt.struc_metric_1_sep('/data/hcp_test_retest_pp/derivatives/chap', 100)
np.save('/data/hcp_test_retest_pp/derivatives/chap_analysis/struc_metrics/struc_metric_1_sep', sep_forreal)
#show plot--they don't separate! (at least not before 100 evs)

#fisher transform of pearson r then t-test
#

'''patrick sparsity metric'''
#for three subjects: within subj. = ~30; across subjs = ~45
#takes a while to run but code is there

metric_2_forreal = t_rt.struc_metric_2('/data/hcp_test_retest_pp/derivatives/chap')
within = np.array(metric_2_forreal[1])
np.save('/data/hcp_test_retest_pp/derivatives/chap_analysis/struc_metrics/struc_metric_2_within', within)


'''patrick second sparsity metric
for one subject: within subj = 1; across subj. = 1...
'''
v1 = np.delete(np.load(f'{chap_out}/sub-103818/ses-test/vecs.npy'), 0, axis=1)
v2 = np.delete(np.load(f'{chap_out}/sub-103818/ses-retest/vecs.npy'), 0, axis=1)
dcp.subspace_distance_eff(v1,v2)

'''spectra'''
chap_out = '/Users/bwinston/Downloads/chap_out_test'
#chap_out = '/data/hcp_test_retest_pp/derivatives/chap'
global s
s = {}
for sub in ['111312', '105923', '103818']:
    s[sub] = {}
    for ses in ['test','retest']:
        s[sub][ses] = {}
        s[sub][ses]['power'] = {}
        s[sub][ses]['energy'] = {}
        s[sub][ses]['recon'] = np.load(f'{chap_out}/sub-{sub}/ses-{ses}/func/REST1/reconspectra/sub-{sub}_ses-{ses}_task-rest1_dynamic_reconstruction_spectrum.npy')
        s[sub][ses]['recon'] = np.delete(s[sub][ses]['recon'], 0, axis=0)
        for spec in ['power', 'energy']:
            for t in ['mean', 'dynamic']:
                s[sub][ses][spec][t] = np.load(f'{chap_out}/sub-{sub}/ses-{ses}/func/REST1/{spec}spectra/sub-{sub}_ses-{ses}_task-rest1_{t}_{spec}_spectrum.npy')
                s[sub][ses][spec][t] = np.delete(s[sub][ses][spec][t], 0, axis=0)
        s[sub][ses]['power']['normalized'] = np.load(f'{chap_out}/sub-{sub}/ses-{ses}/func/REST1/powerspectra/sub-{sub}_ses-{ses}_task-rest1_normalized_power_spectrum.npy')
        s[sub][ses]['power']['normalized'] = np.delete(s[sub][ses]['power']['normalized'], 0, axis=0)
#unordered correlation btwn spectra within/across subject    
print('within subj. mean energy correlation:' + str(pearsonr(s['103818']['test']['energy']['mean'], s['103818']['retest']['energy']['mean'])[0]))
print('within subj. mean pwr correlation:' + str(pearsonr(s['103818']['test']['power']['mean'], s['103818']['retest']['power']['mean'])[0]))
print('across subj. mean energy correlation:' + str(pearsonr(s['103818']['test']['energy']['mean'], s['105923']['test']['energy']['mean'])[0]))
print('across subj. mean pwr correlation:' + str(pearsonr(s['103818']['test']['power']['mean'], s['105923']['test']['power']['mean'])[0]))

#step 1: optimally order sub1ses1 with sub1ses2
test_retest_rel_2v(f'{chap_out}/sub-103818/ses-test/vecs.npy',f'{chap_out}/sub-103818/ses-retest/vecs.npy',100)
#test_retest_rel_2v(f'{chap_out}/sub-103818/ses-test/vecs.npy',f'{chap_out}/sub-105923/ses-test/vecs.npy',100)
cd['ret_ind'] = []
for ev in range(99):
    cd['ret_ind'].append(cd['pairs'][ev]['ret_ind'])
#ordered correlation btwn spectra within/across subject
print('ordered within subj. mean energy correlation:' + str(pearsonr(s['103818']['test']['energy']['mean'], s['103818']['retest']['energy']['mean'][cd['ret_ind']])[0]))
print('ordered within subj. mean pwr correlation:' + str(pearsonr(s['103818']['test']['power']['mean'], s['103818']['retest']['power']['mean'][cd['ret_ind']])[0]))

print('ordered across subj. mean energy correlation:' + str(pearsonr(s['103818']['test']['energy']['mean'], s['105923']['test']['energy']['mean'][cd['ret_ind']])[0]))
print('ordered across subj. mean pwr correlation:' + str(pearsonr(s['103818']['test']['power']['mean'], s['105923']['test']['power']['mean'][cd['ret_ind']])[0]))

'''
recon spectrum:
histogram of values/frequency of plot (zero crossings)
shannon entropy for histogram (how compressible)

compare REST1 and REST2 within session and also REST1 to REST1 across sessions
'''
#how many zero crossings per harmonic on avg. (plot)
t_rt.ev_freq('/Users/bwinston/Downloads/chap_out_test', 100)

#comparing within vs. across subject for spectrum of frequencies
print(t_rt.freq_comp('/Users/bwinston/Downloads/chap_out_test', 100))
#0.42 to .346

#AM I DOING SHANNON ENTROPY RIGHT? discretize, then run function
t_rt.ev_shan_ent('/Users/bwinston/Downloads/chap_out_test', 100)
#shannon entropy of each ev timeseries. list of 99 entropies, compare that across sessions to get within session, e.g.
#.626 within, .56 across


'''RSN STUFF'''



'''

def test_retest_rel_2v(vec_1, vec_2, n_evecs):
    n_evecs = n_evecs-1
    global cd 
    cd = {} #for chap_data
    cd['bcorrs'] = []
    cd['vec_1'] = np.load(vec_1)
    cd['vec_2'] = np.load(vec_2)
    for i in ['1', '2']:   
        cd[f'vec_{i}'] = np.delete(cd[f'vec_{i}'], 0, axis=1)
    global hp
    hp, hp['holy'] = {}, {} #dict for within these fxns
    hp['corr_orig'] = np.empty((n_evecs,n_evecs))
    for evec_1 in range(0,n_evecs): 
        for evec_2 in range(0,n_evecs): #n_evecs x n_evecs correlation matrix
            hp['corr_orig'][evec_2,evec_1] = abs(pearsonr(cd[f'vec_1'][:,evec_1], cd[f'vec_2'][:,evec_2])[0]) #comparing column with column (ev with ev)
    hp['corr_all'] = hp['corr_orig'].swapaxes(0,1) #prepare to turn into dicts
    hp['corr_all'] = {index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(hp['corr_all'])} #turn into dicts
    find_bcorrs_2v(hp, 0, n_evecs) #run find bcorrs function, get best pairs
    for ev in range(n_evecs):
        cd['bcorrs'].append(cd['pairs'][ev]['bcorr']) #all ideal corrs (w/ no repeats)
    return cd['bcorrs']

def find_bcorrs_2v(hp, run, n_evecs): 
    while len(hp['corr_all']) > 0:
        hp[run] = {}
        hp[run]['maxes'] = {}
        for ev in hp['corr_all']: #for each test session evec, which retest session evec is the best match?
            hp[run]['maxes'][ev] = {} #test session ev's dict
            hp[run]['maxes'][ev]['bcorr'] = max(hp['corr_all'][ev].values()) #best correlation w/ any retest ev
            hp[run]['maxes'][ev]['ret_ind'] = get_key(hp['corr_all'][ev], max(hp['corr_all'][ev].values())) #which retest ev was above?
        hp[run]['ret_used'], hp[run]['win_ind'], hp[run]['good_boys'] = [], [], []
        for ev in hp['corr_all']:
            hp[run]['ret_used'].append(hp[run]['maxes'][ev]['ret_ind']) #tally retest indices used above
        ret_dups = set([x for x in hp[run]['ret_used'] if hp[run]['ret_used'].count(x) > 1]) #retest evs that were used multiple times
        if ret_dups == 0:
            break
        hp[run]['ret_used'] = set(hp[run]['ret_used']) #take out duplicates from ret_used
        for ret_dup in ret_dups: #for each retest ev in demand by multiple
            competition = [] #competition is within the retest ev
            for ev in hp[run]['maxes']: 
                if ret_dup == hp[run]['maxes'][ev]['ret_ind']: #if a test evec is involved in a dispute
                    competition.append(hp[run]['maxes'][ev]['bcorr']) #add its correlation to retest ev's competition 
            competition = max(competition) 
            for ev in hp[run]['maxes']:
                if hp[run]['maxes'][ev]['bcorr'] == competition: #if its correlation was a winner, leave as is
                    hp[run]['win_ind'].append(ev) #add test ev that won to win_ind
        for ev in hp['corr_all']:
            if hp[run]['maxes'][ev]['ret_ind'] not in ret_dups: #if test ev not involved in any disputes, add to good_boys
               hp[run]['good_boys'].append(ev) 
        for ev in list(hp['corr_all']):
            if ev in (set(hp[run]['good_boys']) | set(hp[run]['win_ind'])):
               hp['holy'][ev] = hp[run]['maxes'][ev] #add good boys and winners to holy
               del hp['corr_all'][ev] #delete them from corr_all
        for ev in hp['corr_all']:
            for r_ev in hp[run]['ret_used']:
                del hp['corr_all'][ev][r_ev] #idk what this is doing
        run = run + 1
        find_bcorrs_2v(hp, run, len(hp['corr_all'])) #rerun function on smaller corr_all (leftovers)
    cd['pairs'] = hp['holy']     
'''