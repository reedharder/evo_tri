# -*- coding: utf-8 -*-
"""
Created on Sat Feb 21 15:15:56 2015

@author: Reed
"""
import numpy as np
import pandas as pd
from itertools import combinations
from math import sqrt


import os

os.chdir("C:/Users/Reed/Desktop/minjun_tajima")
#function takes rs number, chromosome, population, window, and hapmap phase, return's Tajima's D
def tajima_calc(rs=12083193, chrom=1, pop="ASW", window=10000, phase="3.3"):
    #get file name
    filename = "genotypes_chr%s_%s_phase%s_consensus.b36_fwd" % (chrom, pop, phase)
    #read file into dataframe
    df = pd.read_csv(filename, delim_whitespace=True)    
    
    #get start and end positions of sequence based on window
    rs_position =int(df[df["rs#"]==("rs" + str(rs))]['pos'])# position of requested
    sequence_bound_lower = rs_position - window
    sequence_bound_upper = rs_position + window
    
    #get all snps within window of sequence, and get columns of resulting data frame that are samples (12th column and onward)
    window_df=df[(df["pos"] >= sequence_bound_lower) & (df["pos"] <= sequence_bound_upper)].iloc[:,11:]
    
    #get number of segregating sites in window
    num_snps = window_df.shape[0] #total number of snps in window
    num_nonsegregating= np.count_nonzero(window_df.eq(window_df.iloc[:, 0], axis=0).all(1))#count all snps where all samples have same genotype, i.e. number of nonsegregating sites 
    S = num_snps - num_nonsegregating #number of segregating sites
    
    #get column index pairs of all combinations of sample pairs    
    n = window_df.shape[1] #number of samples 
    sample_combos = combinations(range(0,n),2) # all pairs of combos
    
    #loop through all sample combinations, count number of polymorphisms in each
    polymorphisms = [] #initialize list of number of polymorphisms for each pair
    for combo in sample_combos: #for each sample pair...
        #test sample pair for equality at all sites, return vector where true if genotype matches, false otherwise
        pair_genotype_match = window_df.iloc[:,combo[0]].eq(window_df.iloc[:,combo[1]])
        #count number of falses for this pair, add to list
        polymorphisms.append(num_snps - np.count_nonzero(pair_genotype_match))
    #calculate average number of polymorphisms over all sample pairs
    k_hat = sum(polymorphisms)/len(polymorphisms)
    #calculate other parameters
    a1 = sum([1/i for i in range(1,n)])
    a2 = sum([1/i**2 for i in range(1,n)])
    b1 = (n + 1)/(3*(n-1))
    b2 = (2*(n**2 + n + 3))/(9*n*(n-1))
    c1 = b1 - 1/a1
    c2 = b2 - (n+2)/(a1*n) + a2/(a1**2)
    e1 = c1/a1
    e2 = c2/(a1**2 + a2)
    #calculate difference between genetic diversity estimates
    d = k_hat -  S/a1 #second term is theta, or M, equal to 4*N*mu
    #calculate standard deviation of d
    stdev_d = sqrt(e1*S + e2*S*(S-1))
    #calculate Tajima's D
    D = d/stdev_d
    
    return D
    

