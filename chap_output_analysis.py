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

ties = set([ind for ind in test_maxes['ind'] if test_maxes['ind'].count(ind) > 1])
ind = np.array(test_maxes['ind'])
dups = []
for tie in ties:
    dups.append(np.where(test_maxes['ind'] == tie)[0])
x = []
for i in range(len(dups)):
    x.append(list(dups[i]))   
troublemakers = []
for sublist in x:
    for item in sublist:
        troublemakers.append(item)

winners = []
for pair in x:
    corrs = []
    for h in pair:
        print('test evec: ' + str(h))
        print('its correlation is: ' + str(test_maxes['corr'][h]))
        corrs.append(test_maxes['corr'][h])
    print('retest evec theyre fighting for is: ' + str(test_maxes['ind'][h]))
    print(corrs)
    winners.append(max(corrs))
    print('')

copy = test_maxes['corr'].copy()
for element in range(len(copy)):
    if element in troublemakers:
        copy[element] = 0

for corr in winners:
    print(test_maxes['corr'].index(corr))
    troublemakers.remove(test_maxes['corr'].index(corr))
    copy[test_maxes['corr'].index(corr)] = test_maxes['corr'][test_maxes['corr'].index(corr)]

used_retest = list(dict.fromkeys(test_maxes['ind']))
leftovers = list(set(range(200))- set(used_retest))

def get_key(my_dict, val):
    for key, value in my_dict.items():
         if val == value:
             return key

'''
important shit
'''

test_retest_rel('/Users/bwinston/Downloads/chap_out_test', 200)

def test_retest_rel(chap_dir, n_evecs):
    global cd 
    cd = {} #for chap_data
    subject_dirs = glob(os.path.join(chap_dir, "sub-*"))
    subs = [subject_dir.split("-")[-1] for subject_dir in subject_dirs] 
    for sub in subs:
        cd[sub] = {}
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

     
#retest evecs have to be dict also       
for ev in hp['corr_all']:
   hp['corr_all'][ev][ev] = {} 
    
    run = run + 1
    
    find_bcorrs(hp, run, corr_mat = )



        #add the test vector whose bcorr corresponds to max(competition)


      
    
    hp[run]['ret_ties'] = set([ind for ind in hp[run]['maxes']['ind'] if hp[run]['maxes']['ind'].count(ind) > 1]) #retest evecs that were chosen twice
    for rt_ev in hp[run]['ind']:
 
    
    
    hp[run]['ind'] = np.array(hp[run]['maxes']['ind']) #list of retest evecs chosen
    hp[run]['dups'] = []
    for tie in hp[run]['ties']:
        hp[run]['dups'].append(np.where(hp[run]['maxes']['ind'] == tie)[0])
    hp[run]['dup_pairs'] = []
    for i in range(len(hp[run]['dups'])):
        hp[run]['dup_pairs'].append(list(hp[run]['dups'][i]))   
    hp[run]['troublemakers'], hp[run]['winners'], hp[run]['win_ind'] = [], [], []
    for pair in hp[run]['dup_pairs']: #formerly x. each pair is list of two test evecs that are fighting over same retest evec
        for ind in pair:
            hp[run]['troublemakers'].append(ind) #flat list of every test evec involved in dup_pairs
            
    for pair in dup_pairs:
        competition = []
        for h in pair: #each of the competing test_evecs
            competition.append(hp[run]['maxes']['bcorr'][h]) #append their correlation to list
        hp[run]['winners'].append(max(competition)) #take highest correlation in pair
    #add to_holy (non-troublemakers) and competition winners to holy:
    hp[run]['to_holy'] = list(set(range(n_evecs)) - set(hp[run]['troublemakers'])) #list of test evecs good to go   
    for good in hp[run]['to_holy']: #good test evec
        hp['holy'][good] = {}
        hp['holy'][good]['corr'] = hp[run]['maxes']['bcorr'][good] #corr of that test evec
        hp['holy'][good]['retest_ind'] = hp[run]['maxes']['ind'][good] #retest ind of that test evec
    for w_corr in hp[run]['winners']:
        hp[run]['win_ind'].append(hp[run]['maxes']['bcorr'].index(w_corr)) #winning test evec indices
    for good in hp[run]['win_ind']:
        hp['holy'][good] = {}
        hp['holy'][good]['corr'] = hp[run]['maxes']['bcorr'][good]
        hp['holy'][good]['retest_ind'] = hp[run]['maxes']['ind'][good]
        
    leftover_mat(hp, n_evecs, run)
    
#make a new array with all leftovers and their evecs. then repeat the process

def leftover_mat(hp, n_evecs, run):
#get leftovers from holy, if no leftovers, exit while loop, if still leftovers, prep to do over again
    hp[run]['moving_on'] = list(set(hp[run]['win_ind']) | set(hp[run]['to_holy']))
    hp[run]['troublemakers'] = list(set(range(n_evecs)) - (set(hp[run]['win_ind']) | set(hp[run]['to_holy'])))
    hp[run]['troublemakers'].sort()
    used_retest = list(dict.fromkeys(hp[run]['maxes']['ind']))
    hp[run]['leftovers'] = list(set(range(n_evecs)) - set(used_retest))
    hp[run]['leftovers'].sort()
    n_evecs = len(hp[run]['troublemakers'])
    if n_evecs == 0:
        #break
        print('hi')
    else:
        hp['corr_temp'] = np.delete(hp['corr_temp'], hp[run]['moving_on'], axis = 0)
        hp['corr_temp'] = np.delete(hp['corr_temp'], used_retest, axis = 1)
        
       
        for test_ev in range(len(hp['corr_temp'])):
           if test_ev not in hp[run]['troublemakers']:
               hp['corr_temp'] = np.delete(hp['corr_temp'], test_ev, axis = 0)
        run = run + 1
        find_bcorrs(hp, run, hp['corr_temp'], n_evecs)  
        
       
           

hi = {index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(hp['corr_all'])}

asdf = list(hi[0].values())

hp['corr_temp'] = np.delete(hp['corr_temp'], 3, axis = 0)
           
find_bcorrs(hp, run, corr_mat, n_evecs = 200):           
       
max(hp['corr_all'][0].values())
    

n_evecs = 200  
global hp
hp = {}
hp['holy'] = {}
hp['vecs_test'] = np.load('/Users/bwinston/Downloads/chap_out_test/sub-200109/ses-test/vecs.npy')
hp['vecs_retest'] = np.load('/Users/bwinston/Downloads/chap_out_test/sub-200109/ses-retest/vecs.npy')
hp['corr_orig'] = np.empty((n_evecs,n_evecs))
for evec_test in range(0,n_evecs): 
    for evec_retest in range(0,n_evecs): #n_evecs x n_evecs correlation matrix
        hp['corr_orig'][evec_retest,evec_test] = abs(pearsonr(hp['vecs_test'][:,evec_test], hp['vecs_retest'][:,evec_retest])[0])
hp['corr_all'] = hp['corr_orig'].swapaxes(0,1)
hp['corr_all'] = {index:{i:j for i,j in enumerate(k) if j} for index,k in enumerate(hp['corr_all'])}
find_bcorrs(200109, hp, 0, hp['corr_all'], n_evecs) 


 
#TODO: add stuff to holy, make the funcitons work in concert and iterated 
    
run = 0
while hp['troublemakers'] > 1:
    find_bcorrs(sub, vecs_test, vecs_retest, n_evecs) 
    break_ties(hp)
    leftover_mat(hp, 200)
    run = run + 1


         

retest_maxes = {}
retest_maxes['corr'],  retest_maxes['ind'] = [], [] #for each retest session evec, which test session evec is the best match?
for corr in range(len(corr_200109[0])):
    retest_maxes['corr'].append(max(corr_200109[corr])) #maximum correlation of that evec to another evec
    retest_maxes['ind'].append(np.where(corr_200109[corr]==retest_maxes['corr'][corr])[0][0]) #which evec from other session does the max correspond to

mean_of_good_fits = []
for x in range(1,200):
    mean_of_good_fits.append(stats.mean(retest_maxes['corr'][:x]))
plt.plot(mean_of_good_fits, 'g')
