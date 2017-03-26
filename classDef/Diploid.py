#!/usr/bin/python
# -*- coding: UTF-8 -*-
from .CrossSite import CrossSite

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
        return repr([self.seq1,self.seq2])
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