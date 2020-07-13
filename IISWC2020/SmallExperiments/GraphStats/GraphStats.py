# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 17:44:59 2020

@author: reetb
"""
import networkx as nx
import os
os.chdir('C:/Local/Coursework/Research-HPCBIO-Lab/GraphOrdering/HiPC2020/SmallExperiments/GraphStats/')

filename = 'cti.edgelist'

directed = False

#if (sys.argv[2] == 1):
#    directed = False

if directed:
    G = nx.read_edgelist(filename, nodetype = int, create_using = nx.DiGraph)
else:
    G = nx.read_edgelist(filename, nodetype = int)
    
    
clusCoeff = nx.average_clustering(G)
print('Clustering Coeff: ', clusCoeff)

try:
    diameter = nx.diameter(G)
#    avgPathLen = nx.average_shortest_path_length(G)
#    print('Average Path Length: ', avgPathLen)
    print('Diameter: ', diameter)
except:
    print('Graph not connected')
alltriangles = nx.triangles(G)
triangles = sum(alltriangles.values()) / 3
print('No of triangles: ', triangles)

