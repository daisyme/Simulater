#!/usr/bin/python
# -*- coding: UTF-8 -*-
import random as rd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from .Population import Population
from .Founder import Founder
from .BasicInfo import Generations, Genome

class Experiment:
    '''
    founderNum: the size of founders
    '''
    def __init__(self, founderNum=100, initGeneration= 10, replicate=5, founder=None, SNP=None, mutation=None, output='out.txt'):
        self.output = open(output, 'w+')
        self.rsize = replicate
        if founder is not None:
            self.fsize = founder.size
            self.genome = founder.genome
            self.founder = founder
        else:
            self.fsize = founderNum
            self.genome = Genome()
            self.founder = Founder()
        self.pop = Population(self.fsize, self.rsize, self.genome.cM)
        self.gsize = 0
        self.generations = Generations(self.pop)
        # if it's too big, we can just output them without storing them
        self.run(generation=initGeneration)
    def run(self, migrationRate=0, generation=10, selectionParameter=None):
    #selectionParameter only add the select info from ancestor
        while generation:
            self.output.write(repr(self.gsize))
            self.output.write('\n')
            self.output.write(repr(self.pop))
            self.output.write('\n')
            generation -= 1
            self.gsize += 1
            gameteid = self.pop.selection(selectionParameter,self.founder)
            self.pop.recombineNmate(gameteid)
            self.pop.migrate(migrationRate)
            self.generations.append(self.pop)
            #generations is the place to store all the processing info
    def summary(self, selectionParameter=None, kwd='default'):
        #TODO
        #fst? heterosity?
        Nsnp = len(selectionParameter)
        for replicate in self.pop:
            geno = [0 for i in range(self.fsize)]
            di = 0
            for diploid in replicate:
                for gamete in [diploid.seq1, diploid.seq2]:
                    ig = 0
                    isnp = 0
                    while isnp < Nsnp:
                        if gamete[ig].crossSite.loc > founder.snp[0][selectionParameter[isnp][0]]:
                            csnp = founder.snp[gamete[ig].CrossSite.pid][selectionParameter[isnp][0]]
                            MAF[isnp] += csnp
                            geno[di] += csnp*selectionParameter[isnp][1]
                            isnp += 1
                        ig += 1
                di += 1
        pass
    def draw(self, generation=None, width=0.1, mode=0):
        #LATER: change it into draw snps
        if generation is None:
            generation = self.gsize
        nFig = plt.figure()
        pop = self.generations.generations[generation].pop
        colormap = plt.cm.gist_ncar   #you can change theme if you want
        color = [colormap(i) for i in np.linspace(0, 1, self.fsize)]
        y = 0
        i = 0
        for replicate in pop:
            width = 1/ self.fsize/2
            i += 1
            nAx = nFig.add_subplot(1,self.rsize,i)
            for diploid in replicate:
                for gamete in [diploid.seq1, diploid.seq2]:
                    pLoc = 0
                    y += width
                    for crossSite in gamete:
                        length = (crossSite.loc - pLoc)/self.genome.cM
                        x = pLoc/self.genome.cM
                        p = patches.Rectangle(
                            (x, y), length,width,
                            facecolor=color[crossSite.aid-1]
                        )
                        #print(x, y, length, width)
                        nAx.add_patch(p)
                        pLoc = crossSite.loc
            y = 0
        #nAx.set_xlim([-0.2,6])
        #nAx.set_ylim([-0.1,1.5])
        plt.show()
    def drawSNP(self, generation=None, width=0.1, mode=0, selectionParameter=None):
        #LATER: change it into draw snps
        if generation is None:
            generation = self.gsize
        nFig = plt.figure()
        pop = self.generations.generations[generation].pop
        colormap = plt.cm.gist_ncar   #you can change theme if you want
        founder = self.founder
        Nsnp = len(selectionParameter)
        color = [colormap(i) for i in np.linspace(0, 1, 2)]
        y = 0
        i = 0
        length = 1/Nsnp
        width = 1/self.fsize/2
        for replicate in pop:
            i += 1
            nAx = nFig.add_subplot(1,self.rsize,i)
            for diploid in replicate:
                for gamete in [diploid.seq1, diploid.seq2]:
                    x = 0
                    ig = 0
                    isnp = 0
                    y += width
                    while isnp < Nsnp:
                        #print(ig, isnp, selectionParameter[isnp][0],gamete[ig].loc)
                        if gamete[ig].loc > founder.snp[0][selectionParameter[isnp][0]]:
                            csnp = founder.snp[gamete[ig].pid+1][selectionParameter[isnp][0]]
                            p = patches.Rectangle(
                                ( x, y), length,width,
                                facecolor=color[int(csnp)]
                            )
                            nAx.add_patch(p)
                            #print(x,y,length,width, color[isnp])
                            x += length
                            isnp += 1
                            ig -= 1
                        ig += 1
            y = 0
        plt.show()