# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 17:36:41 2020

@author: reetb
"""

import networkx as nx
import sys


filename = sys.argv[1]
directed = True

if (sys.argv[2] == 1):
    directed = False

if directed:
    G = nx.read_edgelist(filename, nodetype = int, create_using = nx.DiGraph)
else:
    G = nx.read_edgelist(filename, nodetype = int)

# Relabel nodes to start from zero
G = nx.convert_node_labels_to_integers(G, first_label = 1)    

# Change the 2nd argument to just 'filename' below to rewite over the input file
nx.write_edgelist(G, filename, data = False)