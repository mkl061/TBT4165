#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 18:21:25 2023

@author: Jann
"""
#%% IMPORTING LIBRARIES

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pyvis.network import Network
from scipy import stats as st
from project2.s_core import s_core
import copy

#%% TASK ONE
# Read files
filename = "project2/BioGRID_PIN.txt"
df = pd.read_csv(filename, sep='\t')
# \t denotes that the columns are separated by tabs

# Create a network
G_pin = nx.from_pandas_edgelist(df, source='Systematic Name Interactor A', target='Systematic Name Interactor B',edge_attr=True)

# Remove the self loops
selfloops=list(nx.selfloop_edges(G_pin))
G_pin.remove_edges_from(selfloops)

# Network size and node degree
n_nodes = G_pin.number_of_nodes()
n_edges = G_pin.number_of_edges()
av_nd=round(sum([i[1] for i in G_pin.degree()])/n_nodes,2) 
print("\nThe whole PIN network has {} nodes and {} edges \nand an average node degree of {}.".format(n_nodes, n_edges, av_nd))


#%% TASK TWO
tax_id=559292 #Taxonomy ID for yeast 
nodes_yeast=[]

for a,b in G_pin.edges:
    if G_pin.edges[a,b]["Organism Interactor A"]==tax_id and a not in nodes_yeast:
        nodes_yeast.append(a)
    if G_pin.edges[a,b]["Organism Interactor B"]==tax_id and b not in nodes_yeast:
        nodes_yeast.append(b)

G_yeast=G_pin.subgraph(nodes_yeast) # graph consisting only of nodes that stem from yeast
G_yeast=nx.Graph(G_yeast) # since the subgraph will be frozen and can't be modified further

# removing isolates
isos=list(nx.isolates(G_yeast)) # (no isolated nodes)
G_yeast.remove_nodes_from(isos)


# Network size
n_nodes=G_yeast.number_of_nodes()
n_edges=G_yeast.number_of_edges()
av_nd=round(sum([i[1] for i in G_yeast.degree()])/n_nodes,2) 
print("\nThe yeast PIN network has {} nodes and {} edges \nand an average node degree of {}.".format(n_nodes, n_edges, av_nd))

# degree distribution
degree_freq = nx.degree_histogram(G_yeast)   # list of the frequency of each degree value by index
degrees = range(len(degree_freq)) # indices of degree_freq = node degree
# truncated
ra=70
degree_freq_trunc=degree_freq[0:ra]
degrees_trunc=degrees[0:ra]
# Plot figures
plt.figure(figsize=(12, 8))
x_vals = [degrees, degrees_trunc]
y_vals = [degree_freq, degree_freq_trunc]
titles = ["Whole degree distribution", "Degree distribution from degrees 0 to "+str(ra)]
for i in range(len(x_vals)):
    plt.subplot(1, 2, i + 1)
    plt.bar(x_vals[i], y_vals[i], width=0.80, color='b')
    plt.xlabel("Degree")
    plt.ylabel("Frequency")
    plt.title(titles[i], loc='center', fontsize=12, fontweight='bold')

plt.show()

#%% TASK THREE
# Yeast network
# giant component extraction in order to compute descriptors
Gcc = sorted(nx.connected_components(G_yeast), key=len, reverse=True)  # find network components
G_yeast = G_yeast.subgraph(Gcc[0])
# descriptors
clust_coeff=nx.average_clustering(G_yeast)
sh_path=nx.average_shortest_path_length(G_yeast)
diameter=nx.diameter(G_yeast)
desc=[clust_coeff, sh_path, diameter]
print("\nYeast network properties: ")
print("\n".join("{} {:10.4f}".format(x, y) for x, y in zip(['<C>', '<L>', '<D>'], desc)))

#%% BA Network
# parameters
n,m = n_nodes, round(n_edges/n_nodes) # same  #nodes as yeast network, m= #edges to attach from a single node
sim=1 # number of simulated BA networks in order to calculate the mean of sim BA networks
# network descriptors
clust_coeff=[]
sh_path=[]
diameter=[]

# average over 10 graphs
for i in range(sim):
    G_ba=nx.barabasi_albert_graph(n,m)
    clust_coeff.append(nx.average_clustering(G_ba))
    sh_path.append(nx.average_shortest_path_length(G_ba))
    diameter.append(nx.diameter(G_ba))
desc_means = [np.mean(clust_coeff), np.mean(sh_path), np.mean(diameter)] # average over all 5 networks
print("\nBA network properties: ")
print("\n".join("{} {:10.4f}".format(x, y) for x, y in zip(['<C>', '<L>', '<D>'], desc_means)))

#%% ER Network
n,m = n_nodes, n_edges  # same  #nodes as yeast network, m=total edges
clust_coeff=[]
sh_path=[]
diameter=[]

# average over 10 graphs
for i in range(sim):
    G_er=nx.gnm_random_graph(n,m)
    
    # giant component extraction in order to compute descriptors
    Gcc = sorted(nx.connected_components(G_er), key=len, reverse=True)  # find network components
    G_er = G_er.subgraph(Gcc[0])
    
    clust_coeff.append(nx.average_clustering(G_er))
    sh_path.append(nx.average_shortest_path_length(G_er))
    diameter.append(nx.diameter(G_er))
desc_means = [np.mean(clust_coeff), np.mean(sh_path), np.mean(diameter)] # average over all 10 networks
print("\nER network properties: ")
print("\n".join("{} {:10.4f}".format(x, y) for x, y in zip(['<C>', '<L>', '<D>'], desc_means)))

#%% Interpretation
"""
The BA network underestimates the average clustering coefficient as well as the network diameter.
The average shortest path length however is quite accurately estimated.
The same problematic tendency can be seen in the ER network. <C> is underestimated, <L> is 
even more far off than the BA network, whereas <D> is less severly underestimated.
So in general I would argue that the BA and ER are rather bad models for the yeast network.

"""

#%% TASK 4
"""
The k-core analysis is a method to find the "center"/"core" of network. 
The core method recursively "peels" away nodes, the criterion hereby is the node degree.
So first all nodes with degree 1 are removed, then nodes with k=2 and so forth. 
A k-core consists of all nodes with degree >= k. So the 1-core constist of 
nodes with a degree of >=1. So basically all nodes of the network if there are
no isolated nodes.

"""

#%% TASK 5
G_2core=nx.k_core(G_yeast, 2)
# Network size and node degree
n_nodes = G_2core.number_of_nodes()
n_edges = G_2core.number_of_edges()
av_nd=round(sum([i[1] for i in G_2core.degree()])/n_nodes,2) 
print("\nThe 2-core Yeast PIN network has {} nodes and {} edges \nand an average node degree of {}.".format(n_nodes, n_edges, av_nd))
# increase in average node degree makes sense since we peeled away k=1 nodes

# degree distribution
degree_freq_2core = nx.degree_histogram(G_2core)   # list of the frequency of each degree value by index
degrees_2core = range(len(degree_freq_2core)) # indices of degree_freq = node degree
# truncated
ra=70
degree_freq_trunc_2core=degree_freq_2core[0:ra]
degrees_trunc_2core=degrees_2core[0:ra]
# Plot figures
plt.figure(figsize=(12, 8))
x_val = [degrees_2core, degrees_trunc_2core, degrees, degrees_trunc]
y_val = [degree_freq_2core, degree_freq_trunc_2core, degree_freq, degree_freq_trunc]
titles = ["2-core degree distribution", "2-core degree distribution from degrees 0 to "+str(ra), "Original degree distribution",  "Original degree distribution from degrees 0 to "+str(ra)]
for i in range(len(x_val)):
    plt.subplot(2, 2, i + 1)
    plt.bar(x_val[i], y_val[i], width=0.80, color='b')
    plt.xlabel("Degree")
    plt.ylabel("Frequency")
    plt.title(titles[i], loc='center', fontsize=12, fontweight='bold')

plt.show()

"""
It looks like not much has changed. WHY DOES THE LENGTH OF THE DEGREE CHANGE ?? It shouldn't??!
"""


#%% TASK 6
n_nodes_c=G_yeast.number_of_nodes()# number of nodes in original network
n_edges_c=G_yeast.number_of_edges() 
k=0
k_list=[]
n_nodes_list=[]
n_edges_list=[]
while n_nodes_c != 0: # each iteration takes the core with k + 1 from the last iteration until the number of nodes = 0
    k += 1 # updating k
    n_nodes_c=nx.k_core(G_yeast, k).number_of_nodes()  
    n_edges_c=nx.k_core(G_yeast, k).number_of_edges()
      

    n_nodes_c=nx.k_core(G_yeast, k).number_of_nodes() 
    n_edges_c=nx.k_core(G_yeast, k).number_of_edges()

    # for task 7
    k_list.append(k)
    n_nodes_list.append(n_nodes_c)
    n_edges_list.append(n_edges_c)

 
    
    


# The while loop doesn't continue if there are no nodes anymore in the k-core,
# so the last stored k is the k that gave the first core without any nodes 
# therefore, the last core with nodes is the (k-1)-core
innermost_core = nx.k_core(G_yeast,(k-1))

n_nodes = innermost_core.number_of_nodes()
n_edges = innermost_core.number_of_edges()
av_nd=round(sum([i[1] for i in innermost_core.degree()])/n_nodes,2) 
print("\nThe innermost core of the Yeast network is the " + str(k-1)+"-core and has {} nodes and {} edges \nand an average node degree of {}.".format(n_nodes, n_edges, av_nd))

"""
HOW TO CALCULATE THE MAX. NUMBER OF EDGES BETWEEN THESE NODES?
"""


#%% TASK 7
# Plot figures
plt.figure(figsize=(10, 8))
y_vals = [n_nodes_list, n_edges_list]
y_labels = ["# nodes", "# edges"]
titles = ["# nodes vs. k", "# edges vs. k"]
colors = ["b","r"]
for i in range(2):
    plt.subplot(2, 1, i + 1)
    plt.plot(k_list, y_vals[i], str(colors[i])+'o-')
    plt.xlabel("k-core")
    plt.ylabel(y_labels[i])
    plt.xticks(np.arange(min(k_list), max(k_list)+1, 1)) # specify the x-ticks
    plt.title(titles[i], loc='left', fontsize=18, fontweight='bold')
    
#plt.show()


print(n_edges_list)
"""
The nodes seem to decreaes exponentially while the edges seem to decrease rather in a linear fashion.
"""




#%% TASK 8
k_sec=k_list[-3] # last is k with core with no nodes, -2 = innermost, -3 = second innermost
sec_innermost_core=nx.k_core(G_yeast, k_sec)


nt = Network('1000px', '1000px')    # create pyvis network
nt.from_nx(sec_innermost_core)  # convert from networkx to pyvis
nt.show_buttons(filter_=['physics'])
nt.show('nx.html')  #writes network to html file. Visualize using web-browser (e.g. Chrome works, Safari doesn't)


#%% TASK 9
# Loading data
filename = "normalized_expressions.tsv"
exp_df_raw = pd.read_csv(filename, sep='\t')

# Creating a dictionary with gene_id as key and list of gene expression data as value
dicleng=len(exp_df_raw)
exp_dict={exp_df_raw.loc[i][0]: list(exp_df_raw.loc[i][1:]) for i in range(dicleng)}
# apparently there are genes in the expression data that come up more than once

# Pearson and p-value dictionary for each edge
h=list(G_yeast.edges()) #gets a list of all the edges in the Yeast network
attr={} # attribute dictionary which later can be added to the edges as attribute

# st.pearsonr([1,2],[1,2]) # pearson stats

available_exp_genes=[exp_df_raw.loc[i][0] for i in range(len(exp_df_raw))] # gets a list of all genese for which there are expression data available

for i in range(len(h)):
    gene_id_1=str(h[i][0])
    gene_id_2=str(h[i][1])
    
    if gene_id_1 not in available_exp_genes or gene_id_2 not in available_exp_genes:
        p_value=np.nan
        pearson_stats=np.nan
        attr[h[i]]={"p-value": p_value, "pearson-coeff": pearson_stats}
        
    else:
        p_value=float(st.pearsonr(exp_dict[gene_id_1],exp_dict[gene_id_2])[1])
        pearson_stats=float(st.pearsonr(exp_dict[gene_id_1],exp_dict[gene_id_2])[0])
        attr[h[i]]={"p-value": p_value, "pearson-coeff": pearson_stats}


#%%
# Edges fulfilling p-value < 0.0001 and pos. corelation
# add the attr dictionary as attribute to edges to the Yeast network
nx.set_edge_attributes(G_yeast, attr)
# Example: P-value of a certain edge:
#print(G_yeast.edges()[('YAL014C', 'YOR036W')]["p-value"])

prio_edges=[] # a list with edges having a pos. corr. and p-value below 0.0001

for edge in attr:
    if attr[edge]["pearson-coeff"] > 0:
        if attr[edge]["p-value"] < 0.0001:
            prio_edges.append(edge)

print("\nSo {} gene pairs have a positive correlation and a p-value below 0.0001.".format(len(prio_edges)))


#%% TASK 10

"""
No, I think since we are looking at binding interactions, we are not interested in the neg. corr. values.
This is because in order to bind, the proteins have to be expressed together. Otherwise, the interaction doesn't make any sense
since they would never "meet".
"""


#%% TASK 11

# Correlation and p-value matrix
yeast_nodes=list(G_yeast.nodes())

# looking if we have nodes in our Yeast PIN that are  not available in the expression data
notincounter=0
for yn in yeast_nodes:
    if yn not in list(exp_df_raw["Gene ID"]):
        notincounter += 1
print(notincounter)
# apparently we are missing 120 genes from our PIN in the expression data
        

num_ynodes=len(yeast_nodes)
corr_mat= np.empty((num_ynodes,num_ynodes))
corr_mat[:]=np.nan
p_mat= np.empty((num_ynodes,num_ynodes))
p_mat[:]=np.nan


#%%
for column, yncol in enumerate(yeast_nodes):
    for row, ynrow in enumerate(yeast_nodes):
        if yncol in exp_dict.keys() and ynrow in exp_dict.keys():
            corr=float(st.pearsonr(exp_dict[yncol],exp_dict[ynrow])[0])
            p_val=float(st.pearsonr(exp_dict[yncol],exp_dict[ynrow])[1])
            corr_mat[row, column]=corr
            p_mat[row, column]=p_val
            
#%%
G_yeast_s_core=nx.Graph(G_yeast) # again problem with frozen graphs

#%%
# creating a dictionary
# key = s, value=s-core
sl=[2.5,3,3.5,4]
s_core_dic={s:s_core.s_core(G_yeast_s_core, s, yeast_nodes, corr_mat, p_mat) for s in sl }

for s in sl:
    num_edges=s_core_dic[s].number_of_edges()
    num_nodes=s_core_dic[s].number_of_nodes()
    print("\nThe yeast {}-s-core consists of {} edges and {} nodes.".format(s, num_edges, num_nodes))
    #Plot degree
    degree_freq = nx.degree_histogram(s_core_dic[s])   # list of the frequency of each degree value by index
    degrees = range(len(degree_freq)) # indices of degree_freq = node degree
    plt.figure(figsize=(12, 8))
    plt.bar(degrees, degree_freq, width=0.80, color='b')
    plt.xlabel("Degree")
    plt.ylabel("Frequency")
    plt.title("Degree frequency {}-core".format(s), loc='center', fontsize=12, fontweight='bold')


#%% TASK 12
"""
The second inner most core consists of 42 nodes, so the 4-s-core is the best s-core of all the tested
"""
#%% TASK 13
nodes_list=[]
color_list=[]
color=""

for i, node in enumerate(list(G_yeast.nodes)):
    if node in list(s_core_dic[4].nodes):
        color="green"
    if node in list(sec_innermost_core.nodes):
        color="blue"
    if node in list(s_core_dic[4].nodes) and node in list(sec_innermost_core.nodes):
        color="red"
    if color!="":
        color_list.append(color)
        nodes_list.append(node)
    color=""
print(color_list)
    
color_dic={node: color_list[i] for i, node in enumerate(nodes_list)}

hopefully_my_last_network_for_this_exercise=G_yeast.subgraph(nodes_list)

print(hopefully_my_last_network_for_this_exercise)

nx.set_node_attributes(hopefully_my_last_network_for_this_exercise, color_dic, name="color")

print(nx.get_node_attributes(hopefully_my_last_network_for_this_exercise, "color"))


nt = Network('1000px', '1000px')    # create pyvis network
nt.from_nx(hopefully_my_last_network_for_this_exercise)  # convert from networkx to pyvis
nt.show_buttons(filter_=['physics'])
nt.show('last_network.html')  #writes network to html file. Visualize using web-browser (e.g. Chrome works, Safari doesn't)

