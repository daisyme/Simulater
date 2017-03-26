#!/usr/bin/python
# -*- coding: UTF-8 -*-
import random as rd
import numpy as np
from .Diploid import Diploid

class Population:
    '''
        pop is [5][100] and index from 0
        ! put founder info in selectionParameter
    '''
    def __init__(self, founderSize, replicateSize, genomeLength):
        self.fsize = founderSize
        self.rsize = replicateSize
        self.glen = genomeLength
        self.pop = [[Diploid(fi, self.glen) for fi in range(self.fsize)] for ri in range(self.rsize)]
        #how to init sex?
    def __repr__(self):
        return repr(self.pop)
    def recombineNmate(self, gameteid):
        gametePool = [[None for fi in range(self.fsize)] for ri in range(self.rsize)]
        for ri in range(self.rsize):
            for fi in range(self.fsize):
                pid = gameteid[ri][fi]
                crossNum = np.random.poisson(self.glen/100)
                if crossNum == 0:
                    gametePool[ri][fi] = self.pop[ri][pid].recombine(pid)
                while crossNum:
                    cross = rd.uniform(0,self.glen)
                    gametePool[ri][fi] = self.pop[ri][pid].recombine(pid, cross)
                    crossNum -= 1
        #Mate:
        for ri in range(self.rsize):
            rd.shuffle(gametePool[ri])
            #maybe i don't need to shuffle here? since the selection part has already did the shuffle once
            #And didn't pay attention to sex here...sex control is somewhat hard...
            for mi in range(int(self.fsize/2)):
                self.pop[ri][mi*2] = Diploid(seq1=gametePool[ri][mi*2].seq1,seq2=gametePool[ri][mi*2+1].seq2)
                self.pop[ri][mi*2+1] = Diploid(seq1=gametePool[ri][mi*2+1].seq1,seq2=gametePool[ri][mi*2].seq2)        
    def selection(self, selectionParameter=None, founder=None, heritability=0.12, polyHeritability=0.38, optimum=15, std=12):
        '''
        output is a list of selected founders to have the offspring
        my idea is to give one with higher opportunity to have offspring to come up more than once
        '''
        gameteid = []
        if selectionParameter is None:
            for ri in range(self.rsize):
                selected = [i for i in range(self.fsize)]
                rd.shuffle(selected)
                gameteid.append(selected)
        else:
            Nsnp = len(selectionParameter)
            polyMean = 0
            for replicate in self.pop:
                MAF = [0 for i in range(Nsnp)]
                geno = [0 for i in range(self.fsize)]
                di = 0
                for diploid in replicate:
                    for gamete in [diploid.seq1, diploid.seq2]:
                        ig = 0
                        isnp = 0
                        while isnp < Nsnp:
                            #print(ig, isnp, selectionParameter[isnp][0],gamete[ig].loc)
                            if gamete[ig].loc > founder.snp[0][selectionParameter[isnp][0]]:
                                csnp = founder.snp[gamete[ig].pid+1][selectionParameter[isnp][0]]
                                MAF[isnp] += csnp
                                geno[di] += csnp*selectionParameter[isnp][1]
                                isnp += 1
                                ig -= 1
                            ig += 1
                    di += 1
                #?rescale problem
                Vg3 = np.var(geno)
                Ve = np.square(1 - heritability - polyHeritability)
                Vgp = polyHeritability
                pheno = [i for i in range(self.fsize)]
                fitness = [i for i in range(self.fsize)]
                for i in range(self.fsize):
                    pheno[i] = geno[i] + polyMean + rd.gauss(0,np.sqrt(Ve)) + rd.gauss(0,np.sqrt(Vgp))
                    fitness[i] = np.exp(-np.square(pheno[i] - optimum)/np.square(std))
                fsum = np.sum(fitness)
                for i in range(self.fsize):
                    fitness[i] /= fsum
                selected = np.random.choice(self.fsize, self.fsize, p=fitness)
                S = np.mean(np.array(pheno)[selected])
                polyMean += polyHeritability*S
                selected = selected.tolist()
                rd.shuffle(selected)
                gameteid.append(selected)
                print(np.round((np.array(MAF)/self.fsize),2).tolist(), np.mean(pheno), np.sqrt(np.var(pheno)))
            print("---")
        return gameteid                    
    def migrate(self, migrationRate):
        migrationNum = randomMigrate(migrationRate)
        listM = [i for i in range(self.rsize)]
        while(migrationNum):
            migrationNum -= 1
            rd.shuffle(listM)
            self.pop[listM[0]][rd.randint(0,self.fsize-1)] = self.pop[listM[1]][rd.randint(0,self.fsize-1)]
            #!suddenly realize why they record parents relationship and not migration... migration is hard to record relationship between replicates
def randomMigrate(migrationRate):
    return int(abs(np.random.normal(migrationRate,2,1)))
