#!/usr/bin/python
# -*- coding: UTF-8 -*-
from classDef import *
'''
TODO : something haven't done yet
MORE: something want to add but maynot be necessary
?: questions
!: weird thinking
'''
'''
Maybe changing all the data type to array.array for more compact storage?
And rewrite the whole program to pandas.DataFrame can give many convieniece in coding...
'''

F1 = Founder.Founder(rate=0.01)
ExpA = Experiment.Experiment(initGeneration=10,output="out.txt", founder=F1)


selectionParameter2 = [[3,0.1],[5,0.9],[7,0.2],[4,10]]
ExpB = Experiment.Experiment(initGeneration=10,output="out.txt", founder=F1)
ExpB.run(migrationRate=0.2, selectionParameter=selectionParameter2,generation=100)
ExpB.drawSNP(selectionParameter=selectionParameter2)
ExpB.draw()