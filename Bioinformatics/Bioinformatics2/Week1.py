import networkx as nx
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




from collections import defaultdict
'''
Code Challenge: Solve the Overlap Graph Problem (restated below).

Input: A collection Patterns of k-mers.
Output: The overlap graph Overlap(Patterns), in the form of an adjacency list. (You may return the nodes and their edges in any order.)
'''


def overlapgraph(seqs):
	#seqs = list(set(seqs))
	k = len(seqs[0])
	basemap = defaultdict()
	for i, item in enumerate(seqs):
		basemap[item] = []
		for z in seqs[:i] + seqs[i+1:]:
			#print(item[1:])
			if z[:k-1] == item[1:]:
				basemap[item].append(z)

		if basemap[item] == []:
			del basemap[item]


	return(basemap)



'''

Generate a DebrujinGraph Adjacency List representation of a DNA string with k-mer size k
#In Debrujin, Edges represent a k-mer

Inputs: A DNA string and an integer k-mer length k

Output: An adjacency list
'''

def DebrujinGraph(seqs, k):
	G = nx.MultiDiGraph()
	basemap = defaultdict()
	for i in range(len(seqs) - k + 1):
		pref = seqs[i : i + k - 1]
		suf = seqs[i + 1 : i + k]

		G.add_edge(pref, suf, kmer = stringwalk([pref,suf]))
		# if pref in basemap.keys():
		# 	basemap[pref].append(suf)
		# else:
		# 	basemap[pref] = [suf]

	return(G)



  # positions for all nodes - seed for reproducibility

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

G = DebrujinGraph('AAGATTCTCTAAGA', 4, visualize = True)
plt.show()



'''
Generate a DebrujinGraph Adjacency List representation of a list of k-mers

Inputs: A list of k-mers

Output: An adjacency list
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

z = (DebrujinGraph_from_kmer(['GAGG', 'CAGG', 'GGGG', 'GGGA', 'CAGG', 'AGGG', 'GGAG'], visualize = False))
print(z.edges)

plt.show()





