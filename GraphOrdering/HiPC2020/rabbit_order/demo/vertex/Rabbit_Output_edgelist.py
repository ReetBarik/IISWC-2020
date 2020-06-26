# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 13:08:50 2020

@author: reetb
"""
import sys
import networkx as nx

vertexList = sys.argv[1]
inputGraph = sys.argv[2]
directed = True

if (sys.argv[3] == 1):
    directed = False

if directed:
    G = nx.read_edgelist(inputGraph, nodetype = int, create_using = nx.DiGraph)
else:
    G = nx.read_edgelist(inputGraph, nodetype = int)

def read_integers(filename):
    with open(filename) as f:
        return [int(x) for x in f]
    
    
vertices = read_integers(vertexList)

mapping = {}

for i in range(len(vertices)):
    mapping[i] = vertices[i]
    
G = nx.relabel_nodes(G, mapping)

nx.write_edgelist(G, inputGraph, data = False)
