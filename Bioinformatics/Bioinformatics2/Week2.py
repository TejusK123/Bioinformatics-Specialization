import numpy as np 
import networkx as nx 
import matplotlib.pyplot as plt
from collections import deque
import random 
import copy


'''
General Helper Functions:
Import Graph - 
Imports a txt file with a colon seperated adjacency list into python:
Ex. 
0: 3
1: 0
9: 6

gapped_pattern_clean -
Cleans and imports Gapped Pattern Text Input
Ex. 
"GAGA|TTGA TCGT|GATG CGTG|ATGT"

'''
def import_graph(filepath):
	adjacency_list = {}
	with open(filepath, 'r') as f:
		lines = f.readlines()
		
		for line in lines:

			colon_split = line.split(':')
			colon_split[1] = colon_split[1][:-1] if '\n' in colon_split[1] else colon_split[1]
			adjacency_list[colon_split[0]] = deque(colon_split[1][1:].split(' '))

	return(adjacency_list)


def gapped_pattern_clean(GappedPatternpath):
	with open(GappedPatternpath) as f:
		lines = f.readlines()
		k, d = int(lines[0].split(' ')[0]), int(lines[0].split(' ')[1][:-1])
		
		gapped_patterns = [item.split('|') for item in lines[1].split(' ')]
		prefix_kmers = [item[0] for item in gapped_patterns]
		suffix_kmers = [item[1] for item in gapped_patterns]
		suffix_kmers[-1] = suffix_kmers[-1][:-1]
	return(prefix_kmers, suffix_kmers, k, d)




'''
Code Challenge: Solve the Eulerian Cycle Problem.

     Input: The adjacency list of an Eulerian directed graph.
     Output: An Eulerian cycle in this graph.

'''
def Eulerian_Cycle(adjacency_list, solve_type = 'library'):
	if solve_type == 'library':
		G = nx.MultiDiGraph(adjacency_list)
		circuit_edges = list(nx.eulerian_circuit(G))

		return([u for u,v in circuit_edges] + [circuit_edges[-1][-1]])


	start_node = random.choice(list(adjacency_list.keys()))
	start_edges = adjacency_list[start_node]
	stack = deque()
	stack.append((start_node,start_edges))
	
	
	ans = []
	cur_node = start_node
	cur_edges = start_edges
	while len(stack) != 0:
		V = stack[0][1]
		if len(V) == 0:
			ans.append(stack[0][0])
			stack.popleft()
		else:
			cur_node = V.popleft()
			cur_edges = adjacency_list[cur_node]
			stack.appendleft((cur_node,cur_edges))


	return((ans[::-1]))


'''
Helper Function:
num_in_out:
Find Nodes with Unbalanced in-degrees and out-degrees

rotate:
rotates an array 
'''

def num_in_out(adjacency_list):
	outs = [(item,len(adjacency_list[item])) for item in adjacency_list.keys()]
	ins = []
	for item in adjacency_list.values():
		while len(item) != 0:
			ins.append(item.pop())

	
	in_count = [(item,ins.count(item)) for item in set(ins)]
	in_count = sorted(in_count)
	outs = sorted(outs)
	outsdup = outs.copy()
	for item in outs:
		if item in in_count:
			in_count.remove(item)
			outsdup.remove(item)
	


	
	for item in outsdup:
		for item2 in in_count:
			if item[0] == item2[0]:
				in_count.remove(item2)

	try:
		return((in_count[0][0], outsdup[0][0]), True)
	except:
		return((outsdup[0][0], outsdup[-1][0]), False)


def rotate(stuff):
	return(stuff[1:] + stuff[:1])


'''
Code Challenge: Solve the Eulerian Path Problem.

Input: The adjacency list of a directed graph that has an Eulerian path.
Output: An Eulerian path in this graph.

Note: This algorithm is incredibly fucked
'''
def Eulerian_Path(adjacency_list, solve_type = 'library'):
	
	if solve_type == 'library':
		G = nx.MultiDiGraph(adjacency_list)
		path_edges = list(nx.eulerian_path(G))
		return([u for u,v in path_edges] + [path_edges[-1][-1]])
	z = copy.deepcopy(adjacency_list)

	lastedge, dflag = num_in_out(z)
	
	
	flag = False
	
	ins = []
	for item in adjacency_list.values():
		for i in range(len(item)):
			ins.append(item[i])
	
	if dflag == True:
		adjacency_list[lastedge[0]] = deque([lastedge[1]])
		start_node = random.choice(list(adjacency_list.keys()))


	elif ins.count(lastedge[0]) > ins.count(lastedge[1]):
	
		try:
			adjacency_list[lastedge[0]].append(lastedge[1])
			start_node = lastedge[1]
			end_node = lastedge[0]
		except:
			adjacency_list[lastedge[0]] = deque([lastedge[1]])
			start_node = random.choice(list(adjacency_list.keys()))
			flag = True
	else:
		try:
			adjacency_list[lastedge[1]].append(lastedge[0])
			start_node = lastedge[0]
			end_node =lastedge[1]

		except:
			adjacency_list[lastedge[0]] = deque([lastedge[1]])
			start_node = random.choice(list(adjacency_list.keys()))
			flag = True

	if dflag == True:
		start_edges = adjacency_list[start_node]
		stack = deque()
		stack.append((start_node,start_edges))
		
		
		ans = []
		cur_node = start_node
		cur_edges = start_edges
		while len(stack) != 0:
			
			V = stack[0][1]
			
			if len(V) == 0:
				ans.append(stack[0][0])
				stack.popleft()
			else:
				cur_node = V.popleft()
				cur_edges = adjacency_list[cur_node]
				stack.appendleft((cur_node,cur_edges))
			
		ans = ans[::-1][:-1]
	else:
		ans = Eulerian_Cycle(adjacency_list, solve_type = 'na')[:-1]
	if dflag == True:
		while ans[len(ans) - 1] != lastedge[0]:
			ans = rotate(ans)
	else:
		while ans[0] != start_node or ans[-1] != end_node:
			ans = rotate(ans)
	return(ans)






def DebrujinGraph_from_kmer(seqs, visualize = False, solve_type = 'library'):
	if solve_type == 'library':
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
	else:
		basemap = {}
		k = len(seqs[0]) - 1
		for item in seqs:
			pref = item[:k]
			suf = item [len(item) - k:]
			if pref in basemap:
				basemap[pref].append(suf)
			else:
				basemap[pref] = deque([suf])

		return(basemap)
	




def stringwalk(seqs):
	ans = ''
	ans += seqs[0]
	for j in range(1,len(seqs)):
		ans += seqs[j][-1]

	return(ans)

'''
Code Challenge: Solve the String Reconstruction Problem.

Input: An integer k followed by a list of k-mers Patterns.
Output: A string Text with k-mer composition equal to Patterns. (If multiple answers exist, you may return any one.)

'''

def StringReconstruction(Patterns):
	adj_list = DebrujinGraph_from_kmer(Patterns, solve_type = 'manual')
	path = Eulerian_Path(adj_list)
	text = stringwalk(path)
	return(text)





'''
Helper Function:
generates all permuations of length k that can be created from a set of characters
'''

def generate_all_kmers(nucleotide_options, k):

	permutations = [nucleotide for nucleotide in nucleotide_options]

	for i in range(k - 1):
		new_permutations = []
		for permutation in permutations:
			for nucleotide in nucleotide_options:
				new_permutations.append(permutation + nucleotide)
		permutations = new_permutations

	return(permutations)


'''
Code Challenge: Solve the k-Universal Circular String Problem.

Input: An integer k.
Output: A k-universal circular string.
'''

def k_universal_circular_string(k):
	possible_kmers = generate_all_kmers('01', k)
	Debruijin_kmers = DebrujinGraph_from_kmer(possible_kmers, solve_type = 'manual')
	Universal_cycle = Eulerian_Cycle(Debruijin_kmers)
	circular_string = stringwalk(Universal_cycle)
	return(circular_string[:-(k-1)])



'''
Exercise Break: Generate the (3,2)-mer composition of TAATGCCATGGGATGTT in lexicographic order.
Include repeats, and return your answer as a list on a single line.  As a hint to help you with formatting, your answer should begin "(AAT|CAT) (ATG|ATG)..."

'''

def k_d_composition(Seq,k,d):
	ans = sorted([(Seq[i:i+k], Seq[i+d+k:i+d+2*k]) for i in range(len(Seq)-(2*k+d)+1)])
	ans = [f"({item[0]}|{item[1]})" for item in ans]
	return(ans)


'''
Helper Functions:
split_path - 
Splits a path created from a Gapped Pattern Debruijin Graph into prefix and suffix strings

Debrujin_from_Gapped_Pattern:
Creates a Debrujin Graph for Gapped Pattern DNA Strings
'''

def split_path(seq):
	k_meractual_size = int(len(seq[0])/2)
	
	tops = [item[:k_meractual_size] for item in seq]
	bottom = [item[k_meractual_size: ] for item in seq]
	return(tops, bottom)

def Debrujin_from_Gapped_Patten(prefix_kmers, suffix_kmers, Gappedseq = None):
	basemap = {}
	#prefix_kmers, suffix_kmers = ['GAGA', 'TCGT', 'CGTG', 'TGGT', 'GTGA', 'GTGG', 'TGAG', 'GGTC', 'GTCG'], ['TTGA', 'GATG', 'ATGT', 'TGAG', 'TGTT', 'GTGA', 'GTTG', 'GAGA', 'AGAT']
	#prefix_kmers, suffix_kmers = ['AG', 'GC', 'CA', 'AG', 'GC', 'CT', 'TG', 'GC', 'CT'], ['AG', 'GC', 'CT', 'TG', 'GC', 'CT', 'TG', 'GC', 'CA']
	k = len(prefix_kmers[0]) - 1
	for i in range(len(prefix_kmers)):
		item1 = prefix_kmers[i]
		item2 = suffix_kmers[i]

		prefpref = item1[:k]
		sufpref = item1[len(item1) - k:]
		prefsuf = item2[:k]
		sufsuf = item2[len(item2) - k:]

		pref = prefpref + prefsuf
		suf = sufpref + sufsuf
		if pref in basemap:
			basemap[pref].append(suf)
		else:
			basemap[pref] = deque([suf])
	return(basemap)




'''
Code Challenge: Solve the String Reconstruction from Read-Pairs Problem.

Input: Integers k and d followed by a collection of paired k-mers PairedReads.
Output: A string Text with (k, d)-mer composition equal to PairedReads.
'''


def StringSpelledByGappedPattern(GappedPatternpath):
	prefix_kmers, suffix_kmers, k, d = gapped_pattern_clean(GappedPatternpath)
	Debruijin = Debrujin_from_Gapped_Patten(prefix_kmers, suffix_kmers)
	Path = Eulerian_Path(Debruijin)
	prefix, suffix = split_path(Path)
	prefix, suffix = stringwalk(prefix), stringwalk(suffix)
	for i in range(k + d + 1, len(prefix)):
		if prefix[i] != suffix[i - k - d]:
			return("no consensus string can be found")
	
	return(prefix + suffix[-(k+d):])


'''
Helper function:
Generates a list of all nodes with either balanced edges or inbalanced edges

'''


def Balanced_ornot_Degree(adjacency_list, solve_type = 'library'):
	if solve_type == 'library':
		G = nx.MultiDiGraph(adjacency_list)
		one_in_one_out = []
		not_one_in_one_out = []
		for node in G.nodes:
			if len(list(G.successors(node))) != len(list(G.predecessors(node))):
				not_one_in_one_out.append(node)
			else:
				one_in_one_out.append(node)

	else:
		outs = [(item,len(adjacency_list[item])) for item in adjacency_list.keys()]
		ins = []
		for item in adjacency_list.values():
			while len(item) != 0:
				ins.append(item.pop())
		in_count = [(item,ins.count(item)) for item in set(ins)]

		#print(f"\nnum in: {in_count}, num out: {outs}\n")

		one_in_one_out = []
		for item in in_count:
			if item[1] == 1:
				for item2 in outs:
					if (item2[0] == item[0]) and item2[1] == 1:
						one_in_one_out.append(item[0])
		
		not_one_in_one_out = [item for item in list(adjacency_list.keys()) if item not in one_in_one_out]

	return(one_in_one_out, not_one_in_one_out)



'''
Code Challenge: Implement MaximalNonBranchingPaths.

Input: The adjacency list of a graph whose nodes are integers.
Output: The collection of all maximal nonbranching paths in this graph.

#try changing solve type if the function doesn't work
'''
def MaximalNonBranchingPath(adjacency_list, solve_type = 'library'):
	seqs = adjacency_list
	seqscopy = copy.deepcopy(seqs)
	one_in_one_out, not_one_in_one_out = Balanced_ornot_Degree(seqscopy, solve_type)
	#print(not_one_in_one_out)
	paths = []
	for item in list(seqs.keys()):
		if item in not_one_in_one_out:
			#print(f"not one_in_one_out {item}")
			if len(seqs[item]) > 0:
				for i in range(len(seqs[item])):
					cur_node = item
					edges = seqs[cur_node]
					Non_branching = [cur_node, seqs[cur_node][i]]
					cur_node = seqs[cur_node][i]
					while cur_node in one_in_one_out:
						Non_branching.append(seqs[cur_node][0])
						cur_node = seqs[cur_node][0]

					paths.append(Non_branching)

	cycle_paths = []
	for item in list(seqs.keys()):
		if item in one_in_one_out:
			cur_node = item
			edges = seqs[cur_node]
			cur_cycle = [cur_node, seqs[cur_node][0]]
			cur_node = seqs[cur_node][0]
			while cur_node in one_in_one_out and (cur_cycle[0] != cur_cycle[len(cur_cycle) - 1]):
				cur_cycle.append(seqs[cur_node][0])
				cur_node = seqs[cur_node][0]
			
			if (cur_cycle[0] == cur_cycle[len(cur_cycle) - 1]) and (set(cur_cycle) not in [set(item) for item in cycle_paths]):
				cycle_paths.append(cur_cycle)

	#print(cycle_paths)
	
	return(paths + cycle_paths)





'''
Contig Generation Problem: Generate the contigs from a collection of reads (with imperfect coverage).

Input: A collection of k-mers Patterns.
Output: All contigs in DeBruijn(Patterns).
Code Challenge: Solve the Contig Generation Problem.
'''


def Contig_Generation(k_mers):
	Debruijin = DebrujinGraph_from_kmer(k_mers, solve_type = 'manual')
	Non_branching_paths = MaximalNonBranchingPath(Debruijin)
	contigs = [stringwalk(item) for item in Non_branching_paths]
	return(contigs)
























