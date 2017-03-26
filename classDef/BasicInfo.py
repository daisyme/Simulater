#!/usr/bin/python
# -*- coding: UTF-8 -*-
import random as rd
import numpy as np
from .Population import Population

#Idea cM condition
class Chromosome:
    def __init__(self, chrom, cMstart, cMlength):
      self.str = chrom
      self.cM = cMstart
      self.cMlen = cMlength

class Genome:
    '''
    chromosomes: [str, site] [str, site]
    '''
    def __init__(self, cM=80, genomeType=0, chromosomeNum=0, chromosomes=None, **kwds):
        self.cM = cM
        self.chr = chromosomeNum
        self.type = genomeType
        #if we can extend it by built-in types
        if 'bp' in kwds:
            raise NotImplementedError("bp keyword of genome "
                                      "is not implemented")

class Position:
    def __init__(self, cM=None, bp=None, chrom=None, genome=None):
        if genome is not None:
            self.genome = genome
        else:
            self.genome = Genome()
        if chrom is None:
            if cM is not None:
                self.cM = cM
            elif bp is not None:
                self.bp = bp
                #self.cM = bp2cM(bp)
                #TODO
                raise NotImplementedError("bp transfer "
                                      "is not implemented")
            else:
                raise Error("no input position info?")
        else:
            if cM is not None:
                self.cM = chrom.cM + cM
                #well, not sure how you usually present the distance corresponding to the specific chromosome?
        self.bp = bp


def read_genome(io):
#should be derived from pandas.ExcelFile
    return Genome()


class Generations:
    #TODONEXT
    def __init__(self, population):
        self.generations = []
        self.generations.append(population)
    def append(self, population):
        self.generations.append(population)
    def __repr__(self):
        return repr(self.generations)