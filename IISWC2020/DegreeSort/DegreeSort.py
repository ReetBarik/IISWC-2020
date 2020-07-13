# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:05:40 2020

@author: reetb
"""

import sys
import networkx as nx

filename = sys.argv[1]
directed = True

if (sys.argv[2] == 1):
    directed = False

if directed:
    G = nx.read_edgelist(filename, nodetype = int, create_using = nx.DiGraph)
else:
    G = nx.read_edgelist(filename, nodetype = int)

    
l = sorted(G.degree, key = lambda x: x[1], reverse = True)

mapping = {}

for i in range(len(l)):
    mapping[l[i][0]] = i
    
G = nx.relabel_nodes(G, mapping)

nx.write_edgelist(G, filename, data = False)

