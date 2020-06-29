# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 18:33:44 2020

@author: reetb
"""
import random
import sys
import networkx as nx


inputGraph = sys.argv[1]
directed = True

if (sys.argv[2] == 1):
    directed = False

if directed:
    G = nx.read_edgelist(inputGraph, nodetype = int, create_using = nx.DiGraph)
else:
    G = nx.read_edgelist(inputGraph, nodetype = int)


vertices = list(range(G.order()))

random.shuffle(vertices)

mapping = {}

for i in range(len(vertices)):
    mapping[i] = vertices[i]
    
G = nx.relabel_nodes(G, mapping)

nx.write_edgelist(G, inputGraph, data = False)