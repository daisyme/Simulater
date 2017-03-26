#!/usr/bin/python
# -*- coding: UTF-8 -*-
import random as rd
import numpy as np
from .Population import Population
from .BasicInfo import Genome

class Founder(object):
    '''
    snp: sites info + detailed size*snpsize info for each site
    ! I assume that we get the np.ndarray format "snp" data, otherwise, try using np.frombuffer(snp) ,np.fromiter(snp, dtype=np.int) , np.array(snp)   
    '''
    def __init__(self, size=100, snpsize=10, rate=0.1, snp = None, genome=None):
        if genome is None:
            self.genome = Genome()
        if snp is None:
            self.size = size
            self.snpsize = snpsize
            self.snp = np.concatenate(([np.linspace(0,self.genome.cM, snpsize)],  np.random.binomial(1,rate,[size,snpsize])))
        else:
            self.size = snp.shape(1) - 1
            self.snpsize = snp.shape(0)
            self.snp = snp
#MORE
'''
can add a function to adding more snpinfo to current founder
easily done by using 
addsnp = np.random.binomial(1,0,(size,2))
addloc = [2,3]
add = np.concatenate((addloc,addsnp))
snp = np.concatenate((snp,add), axis=1)
'''