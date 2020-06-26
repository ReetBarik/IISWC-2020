# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 18:22:03 2020

@author: reetb
"""

import os
os.chdir('C:/Local/Coursework/Research-HPCBIO-Lab/Graph Ordering/HiPC2020/')

import numpy as np

filename = 'CommunityDetectionOutput.txt'

communities = np.genfromtxt(filename, dtype = int)

community = {}
sizes = []

for i in range(np.max(communities) + 1):
    l = [j for j, x in enumerate(communities) if x == i ]
    sizes.append(len(l))
    community[i] = l
    
