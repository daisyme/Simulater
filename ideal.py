#!/usr/bin/python
# -*- coding: UTF-8 -*-
import pandas as pd
import random as rd
import copy as cp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
'''
TODO: somthing havn't done yet
?: questions
!: weird thinking
'''
'''
Maybe changing all the data type to array.array for more compact storage?
And rewrite the whole program to pandas.DataFrame can give many convieniece in coding...
'''

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

class Population:
    '''
        pop is [5][100] and index from 0
        ! put founder info in selectionParameter
    '''
    def __init__(self, founderSize, replicateSize, genomeLength):
        self.fsize = founderSize
        self.rsize = replicateSize
        self.glen = genomeLength
        self.pop = [founder.pop for ri in range(replicateSize)]
            #how to init sex?
    def __repr__(self):
        return repr(self.pop)
    def recombineNmate(self, gameteid):
        gametePool = [[None for fi in range(self.fsize)] for ri in range(self.rsize)]
        for ri in range(self.rsize):
            for fi in range(self.fsize):
                cross = randomCross()
                pid = gameteid[ri][fi]
                if cross >= self.glen:
                    gametePool[ri][fi] = self.pop[ri][pid].recombine(pid)
                else:
                    #print("pop")
                    #print(repr(self.pop[ri][pid]))
                    gametePool[ri][fi] = self.pop[ri][pid].recombine(pid, cross)
        #Mate:
        for ri in range(self.rsize):
            rd.shuffle(gametePool[ri])
            #maybe i don't need to shuffle here? since the selection part has already did the shuffle once
            #And didn't pay attention to sex here...sex control is somewhat hard...
            for mi in range(int(self.fsize/2)):
                self.pop[ri][mi*2] = Diploid(seq1=gametePool[ri][mi*2].seq1,seq2=gametePool[ri][mi*2+1].seq2)
                self.pop[ri][mi*2+1] = Diploid(seq1=gametePool[ri][mi*2+1].seq1,seq2=gametePool[ri][mi*2].seq2)        
    def selection(self, selectionParameter=None):
        if selectionParameter is None:
            gameteid = []
            for ri in range(self.rsize):
                selected = [i for i in range(self.fsize)]
                rd.shuffle(selected)
                gameteid.append(selected)
            return gameteid
        else:

        #TODO
    def migrate(self, migrationRate):
        migrationNum = randomMigrate(migrationRate)
        listM = [i for i in range(self.rsize)]
        while(migrationNum):
            migrationNum -= 1
            rd.shuffle(listM)
            self.pop[listM[0]][rd.randint(0,self.fsize-1)] = self.pop[listM[1]][rd.randint(0,self.fsize-1)]
            #!suddenly realize why they record parents relationship and not migration... migration is hard to record relationship between replicates

class Generations:
    #TODO
    def __init__(self, population):
        self.generations = []
        self.generations.append(population)
    def append(self, population):
        self.generations.append(population)
    def __repr__(self):
        return repr(self.generations)
class Experiment:
    '''
    founderNum: the size of founders
    '''
    def __init__(self, founderNum=100, initGeneration= 10, replicate=5, founder=None, SNP=None, mutation=None, output='out.txt'):
        self.output = open(output, 'w+')
        if founder is not None:
            self.fsize = founder.size
            self.pop = Population(self.fsize, self.rsize, self.genome.cM)
            self.genome = founder.genome
        else:
            self.fsize = founderNum
            self.pop = Population(self.fsize, self.rsize, self.genome.cM)
            self.genome = Genome()
        self.gsize = 0
        self.rsize = replicate
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
            gameteid = self.pop.selection(selectionParameter)
            self.pop.recombineNmate(gameteid)
            self.pop.migrate(migrationRate)
            self.generations.append(self.pop)
            #generations is the place to store all the processing info
    def summary(self, kwd='default'):
        #TODO
        #fst? heterosity?
        pass
    def draw(self, generation=None, width=0.1, mode=0):
        #LATER: change it into draw snps
        if generation is None:
            generation = self.gsize
        nFig = plt.figure()
        pop = self.generations.generations[generation].pop
        nAx = nFig.add_subplot(111, aspect='equal')
        colormap = plt.cm.gist_ncar   #you can change theme if you want
        color = [colormap(i) for i in np.linspace(0, 1, self.fsize)]
        y = 0
        for replicate in pop:
            width = 0.9/ self.fsize
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
                        print(x, y, length, width)
                        nAx.add_patch(p)
                        pLoc = crossSite.loc
            y += 1
        plt.show()
def randomMigrate(migrationRate):
    return int(abs(np.random.normal(migrationRate,2,1)))

def randomCross():
    #TODO
    '''
    a little function who gives the random breakpoints according to the cM size
    ? you used runif(1,0,1)here... Is that correct? 
    ! your unit is 100cM, And should we just ignore double site recombination?
    '''
    return rd.uniform(0,100)
#
class Founder(object):
    '''
    snp: sites info + detailed size*snpsize info for each site
    ! I assume that we get the np.ndarray format "snp" data, otherwise, try using np.frombuffer(snp) ,np.fromiter(snp, dtype=np.int) , np.array(snp)   
    '''
    def __init__(self, size=10, snpsize=10, snp = None, genome=None):
        if genome is None:
            self.genome = Genome()
        if snp is None:
            self.size = size
            self.snpsize = snpsize
            self.snp = np.concatenate(([np.linspace(0,self.genome.cM, snpsize)],  np.random.binomial(1,0.1,[size,snpsize])))
        else:
            self.size = snp.shape(1)
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

class Diploid(object):
    def __init__(self, founder=0, length=1, sex=0, seq1=None, seq2=None):
        self.seq1 = [CrossSite(founder, founder, length)]
        self.seq2 = [CrossSite(founder, founder, length)]
        if seq1 is not None:
            self.seq1 = seq1
        if seq2 is not None:
            self.seq2 = seq2
        self.sex = sex
        #stop thinking about sex now?
    def __repr__(self):
        return repr((self.seq1,self.seq2))
    def sort(self):
        self.seq1 = sorted(self.seq1)
        self.seq2 = sorted(self.seq2)
    def recombine(self, pid, loc=None):
        o_seq1 = []
        o_seq2 = []
        for i in self.seq1:
            o_seq1.append(i.nextgen(pid))
        for i in self.seq2:
            o_seq2.append(i.nextgen(pid))
        if loc is None:
            return Diploid(seq1=o_seq1,seq2=o_seq2)
        else:
            cro = CrossSite(pid,0,loc)
            o_seq1.append(cro)
            o_seq2.append(cro)
            o_seq1.sort()
            o_seq2.sort()
            i1 = o_seq1.index(cro)
            i2 = o_seq2.index(cro)
            cro1 = cro.cut(self.seq1[i1])
            cro2 = cro.cut(self.seq2[i2])
            seq1 = o_seq1[0:i1] + [cro1] + o_seq2[i2+1:]
            seq2 = o_seq2[0:i2] + [cro2] + o_seq1[i1+1:]
            return Diploid(seq1=seq1,seq2=seq2)


class CrossSite(object):
    def __init__(self, parent, ancestor, location):
        self.pid = parent
        self.aid = ancestor
        self.loc = location
    def __repr__(self):
        return repr((self.pid, self.aid, self.loc))
    def __lt__(self, other):
        return self.loc < other.loc
    def __gt__(self, other):
        return self.loc > other.loc
    def __eq__(self, other):
        return self.loc == other.loc
    def __le__(self, other):
        return self.loc <= other.loc
    def __ge__(self, other):
        return self.loc >= other.loc
    def __ne__(self, other):
        return self.loc != other.loc
    def cut(self, other):
        return CrossSite(self.pid, other.aid, self.loc)
    def nextgen(self, pid):
        return CrossSite(pid, self.aid, self.loc)

#now the sorted is not within the structure which may slow down the speed
#TODO: improve the structure



A = Experiment(output="out.txt")
A.draw()

ax.set_xlim
set_ylim
