#functions for calulating pairwise Fst from arrays of population data

import numpy as np

'''
CaiWeirFst2pop: calculates W-C Fst using Cai formulation

    n1 and n2 are sample sizes of pops 1 and 2 respectively
    p1 and p2 are the allele frequencies for pops 1 and 2 respectively
'''

def CaiWeirFst2pop(n1,n2,p1,p2):

    n=n1+n2
    nc=n-((n1**2 + n2**2)/(n))
    p_hat=(n1*p1 + n2*p2)/n #average sample frequency of allele
    MSP=(n1*(p1-p_hat)**2) + (n2*(p2-p_hat)**2)
    MSG=(1/((n1-1)+(n2-1)))*(n1*p1*(1-p1) + n2*p2*(1-p2))
    with np.errstate(invalid='ignore'): #ignore divide by zero runtime warnings (nan is output as Fst value)
            Fst=(MSP-MSG)/(MSP+(nc-1)*MSG)

    return Fst

'''
WrightFst: calculates Wright's Fst, using formulation here: http://www.uwyo.edu/dbmcd/popecol/maylects/fst.html

    p1 and p2 are the allele frequencies for pops 1 and 2 respectively
'''

def WrightFst(p1,p2):

    h1 = 2*p1*(1-p1) #HWE heterozygosity for allele in p1
    h2 = 2*p2*(1-p2) #HWE heterozygosity for allele in p2
    Hs = (h1+h2)/2 #average heterozygosity of subpopulations
    P = (p1+p2)/2 #average allele frequency
    Ht = 2*P*(1-P) #HWE heterozygosity of total population
    with np.errstate(invalid='ignore'): #ignore divide by zero runtime warnings (nan is output as Fst value)
        Fst = (Ht-Hs)/Ht # calculate wright's Fst

    return Fst