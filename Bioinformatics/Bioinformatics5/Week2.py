import numpy as np 
import networkx as nx
from itertools import combinations




'''
Code Challenge: Implement HierarchicalClustering.

Input: An integer n, followed by an n x n distance matrix.
Output: The result of applying HierarchicalClustering to this distance matrix (using Davg), with each newly created cluster listed on each line.
'''

'''
Helper Function
computes average distance between two clusters
'''


def Distance(D, Cluster1, Cluster2):
	
	distances = []
	for i in range(len(Cluster1)):
		for j in range(len(Cluster2)):
			
			distances.append(D[Cluster1[i],Cluster2[j]])

	return(sum(distances)/(len(Cluster1) * len(Cluster2)))





def UPGMA(D,n):
	Clusters = [(i,[i]) for i in range(n)]
	T = nx.DiGraph()
	Cluster_Creation_List = []
	for item in Clusters:
		T.add_node(item[0], age = 0)
	node_iter = n 
	while len(Clusters) != 1:
		sorted_clusters = sorted(combinations(Clusters,2), key = lambda x: Distance(D, x[0][1],x[1][1]))[0]
		closest_distance = Distance(D, sorted_clusters[0][1], sorted_clusters[1][1])
		
		C_new = (node_iter, sorted_clusters[0][1] + sorted_clusters[1][1])
		Cluster_Creation_List.append(list(map(lambda x: x+1, C_new[1])))
		T.add_edge(C_new[0], sorted_clusters[0][0])
		T.add_edge(C_new[0], sorted_clusters[1][0])
		T.add_edge(sorted_clusters[0][0], C_new[0])
		T.add_edge(sorted_clusters[1][0], C_new[0])
		T.nodes[C_new[0]]['age'] = closest_distance/2
		Clusters.remove(sorted_clusters[0])
		Clusters.remove(sorted_clusters[1])
		Clusters.append((node_iter, sorted_clusters[0][1] + sorted_clusters[1][1]))
		node_iter += 1
	
	for edge in T.edges:
		T[edge[0]][edge[1]]['weight'] = abs(T.nodes[edge[0]]['age'] - T.nodes[edge[1]]['age'])
	
	
	return(T, Cluster_Creation_List)

#usage

n = 7

'''
7
0.00 0.74 0.85 0.54 0.83 0.92 0.89
0.74 0.00 1.59 1.35 1.20 1.48 1.55
0.85 1.59 0.00 0.63 1.13 0.69 0.73
0.54 1.35 0.63 0.00 0.66 0.43 0.88
0.83 1.20 1.13 0.66 0.00 0.72 0.55
0.92 1.48 0.69 0.43 0.72 0.00 0.80
0.89 1.55 0.73 0.88 0.55 0.80 0.00
'''

D = np.loadtxt(r"testdata.txt", skiprows = 1)

Graph, Cluster_Formations = UPGMA(D,n)

for item in Cluster_Formations:
	print(' '.join([str(num) for num in item]))