#Import the required libraries

import csv
from operator import itemgetter
import networkx as nx
from networkx.algorithms import community 
#This part of networkx, for community detection, needs to be imported separately.
from Bio.PDB import *
import numpy as np
import pandas as pd  
from optparse import OptionParser
import sys, os
import io
from io import StringIO
from pathlib import Path
import os.path
from prody import *

filename = fetchPDB('3wbm')
atoms = parsePDB('3wbm')

sele = atoms.select("name CA or name P")
ca = atoms.select("name P")

a,b,dist = [],[],[]
#calcDistance(sele,sele)
for i in range(0,len(sele)):
    for j in range(0,len(sele)):
        distance = calcDistance(sele[i],sele[j])
        if (distance < 5) & (distance!=0):
            dist.append(calcDistance(sele[i],sele[j]))
            a.append(sele[i].getCoords())
            
            b.append(sele[j])

node_names = []
coordinates = []
set_1 = {}
for i in sele:
    node_names.append(i.getResnum())
    coordinates.append(i.getCoords())
set_1[0] = node_names
set_1[1] = coordinates

set_1

distance_matrix = buildDistMatrix(sele,sele,format='mat')

def merge(list1, list2): 
      
    merged_list = [(list1[i], list2[i]) for i in range(0, len(list1))] 
    return merged_list 
edge_list = merge(node_names,b)

G = nx.Graph() # Initialize a Graph object
G.add_nodes_from(node_names) # Add nodes to the Graph
G.add_edges_from(edge_list) # Add edges to the Graph
#print(nx.info(G)) # Print information about the Graph

node_1, node_2,edge_len = [],[],[]
for i in range(0,len(distance_matrix[0])):
    for j in range(0,len(distance_matrix[0])):
        node_1.append((set_1[0][i])*1)
        node_2.append(set_1[0][j])
        edge_len.append(distance_matrix[i][j])

#print (node_1,node_2,edge_len)
#Useee the below code if you want to visualize network
'''
import matplotlib.pyplot as plt

G = nx.from_pandas_edgelist(df,'node1','node2', edge_attr='edge')
durations = [i['edge'] for i in dict(G.edges).values()]
labels = [i for i in dict(G.nodes).keys()]
labels = {i:i for i in dict(G.nodes).keys()}

fig, ax = plt.subplots(figsize=(12,5))
pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, ax = ax, labels=True)
nx.draw_networkx_edges(G, pos, width=durations, ax=ax)
_ = nx.draw_networkx_labels(G, pos, labels, ax=ax)

'''
