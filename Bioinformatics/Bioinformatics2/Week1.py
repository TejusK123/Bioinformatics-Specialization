
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
'''
Code Challenge: Solve the String Composition Problem.

Input: An integer k and a string Text.
Output: Compositionk(Text) (the k-mers can be provided in any order).

'''

def String_Composition(Seq, k):
	return([Seq[i : i + k] for i in range(len(Seq) - k + 1)])



'''
Code Challenge: Solve the String Spelled by a Genome Path Problem.
'''

def stringwalk(seqs):
	ans = ''
	ans += seqs[0]
	for j in range(1,len(seqs)):
		ans += seqs[j][-1]

	return(ans)




'''
Code Challenge: Solve the Overlap Graph Problem (restated below).

Input: A collection Patterns of k-mers.
Output: The overlap graph Overlap(Patterns), in the form of an adjacency list. (You may return the nodes and their edges in any order.)
'''


def overlapgraph(seqs, visualize = False):
	G = MultiDiGraph()
	k = len(seqs[0])
	for i, item in enumerate(seqs):
		for z in seqs[:i] + seqs[i+1:]:
			if z[:k-1] == item[1:]:
				G.add_edge(item,z, kmer = z[:k-1])
	
	if visualize == True:
		pos = nx.planar_layout(G)
		
		names = {name: name for name in G.nodes}
		nx.draw_networkx_nodes(G, pos, node_color = 'b', node_size = 250, alpha = 1)
		nx.draw_networkx_labels(G,pos,names,font_size=8,font_color='w')
		ax = plt.gca()
		for e in G.edges:
		    ax.annotate("",
		                xy=pos[e[1]], xycoords='data',
		                xytext=pos[e[0]], textcoords='data',
		                arrowprops=dict(arrowstyle="->", color="0",
		                                shrinkA=10, shrinkB=10,
		                                patchA=None, patchB=None,
		                                connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*e[2])
		                                ),
		                                ),
		                
		                )
		plt.axis('off')
		print("Graph Image Created: Call \'plt.show()\'")


	return(G)



'''
Code Challenge: Solve the De Bruijn Graph from a String Problem.

Input: An integer k and a string Text.
Output: DeBruijnk(Text), in the form of an adjacency list.
'''

def DebrujinGraph(seqs,k, visualize = False):
	G = nx.MultiDiGraph()
	for i in range(len(seqs) - k + 1):
		pref = seqs[i : i + k -1]
		suf = seqs[i + 1: i + k]

		G.add_edge(pref, suf, kmer = stringwalk([pref,suf]))

	if visualize == True:
		pos = nx.planar_layout(G)
		
		names = {name: name for name in G.nodes}
		nx.draw_networkx_nodes(G, pos, node_color = 'b', node_size = 250, alpha = 1)
		nx.draw_networkx_labels(G,pos,names,font_size=8,font_color='w')
		ax = plt.gca()
		for e in G.edges:
		    ax.annotate("",
		                xy=pos[e[1]], xycoords='data',
		                xytext=pos[e[0]], textcoords='data',
		                arrowprops=dict(arrowstyle="->", color="0",
		                                shrinkA=10, shrinkB=10,
		                                patchA=None, patchB=None,
		                                connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*e[2])
		                                ),
		                                ),
		                
		                )
		plt.axis('off')
		print("Graph Image Created: Call \'plt.show()\'")

	return(G)




'''
DeBruijn Graph from k-mers Problem: Construct the de Bruijn graph from a set of k-mers.

Input: A collection of k-mers Patterns.
Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).
'''
def DebrujinGraph_from_kmer(seqs, visualize = False):
	G = nx.MultiDiGraph()
	k = len(seqs[0]) - 1
	for item in seqs:
		pref = item[:k]
		suf = item [len(item) - k:]
		G.add_edge(pref, suf, kmer = stringwalk([pref,suf]))
	if visualize == True:
		pos = nx.planar_layout(G)
		
		names = {name: name for name in G.nodes}
		nx.draw_networkx_nodes(G, pos, node_color = 'b', node_size = 250, alpha = 1)
		nx.draw_networkx_labels(G,pos,names,font_size=(24*1/k),font_color='w')
		ax = plt.gca()
		for e in G.edges:
		    ax.annotate("",
		                xy=pos[e[1]], xycoords='data',
		                xytext=pos[e[0]], textcoords='data',
		                arrowprops=dict(arrowstyle="->", color="0",
		                                shrinkA=10, shrinkB=10,
		                                patchA=None, patchB=None,
		                                connectionstyle="arc3,rad=rrr".replace('rrr',str(0.3*e[2])
		                                ),
		                                ),
		                
		                )
		plt.axis('off')
		print("Graph Image Created: Call \'plt.show()\'")
	return(G)







