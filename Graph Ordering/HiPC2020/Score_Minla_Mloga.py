# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 15:01:47 2020

@author: reetb
"""
import sys
import numpy as np
import math

filename = sys.argv[1]
linear = False

if (sys.argv[2] == '1'):
    linear = True
    
edgelist = np.genfromtxt(filename, delimiter=' ')

differenceList = [int(x) for x in abs(np.diff(edgelist))]

if not linear:
    differenceList = [math.log10(x) for x in differenceList]
    
score = np.sum(differenceList)

print(score)