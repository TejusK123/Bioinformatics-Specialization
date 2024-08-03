import numpy as np 
import networkx as nx 


'''
Standard AminoAcids with their corresponding Masses
'''

AminoAcidMass = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113,
                 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
                 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}


'''
Mini AminoAcids for testing
'''

minipeptides = {'X' : 4, 'Z' : 5}









'''
Code Challenge: Construct the graph of a spectrum.

Input: A space-delimited list of integers Spectrum.
Output: Graph(Spectrum).
'''

def Spectrum_Graph(Spectrum, Amino_Mass_Dict):
	revAminoAcidMass = {item[1] : item[0] for item in Amino_Mass_Dict.items()}

	Normalize_Spectrum = np.array([0] + Spectrum)

	diffs = np.tril(np.subtract.outer(Normalize_Spectrum, Normalize_Spectrum))

	G = nx.MultiDiGraph()

	for item in Normalize_Spectrum:
		if item in revAminoAcidMass:
			G.add_edge(0, item, amino_acid = revAminoAcidMass[item])

	for i in range(diffs.shape[0]):
		for j in range(diffs.shape[0]):
			if diffs[i,j] in revAminoAcidMass.keys():
				G.add_edge(Normalize_Spectrum[j], Normalize_Spectrum[i], amino_acid = revAminoAcidMass[diffs[i,j]])


	return(G)


#usage
G = Spectrum_Graph([57,71,154,185,301,332,415,429,486], AminoAcidMass)
ans = nx.to_dict_of_dicts(G)



for item in ans.items():
    for i in ans[item[0]]:
        print(f"{item[0]}->{i}: {list((ans[item[0]][i][0].values()))[0]}")






'''
Code Challenge: Solve the Decoding an Ideal Spectrum Problem.

Input: A space-delimited list of integers Spectrum.
Output: An amino acid string that explains Spectrum.
'''



def DecodingIdealSpectrum(Spectrum, Amino_Mass_Dict):

	G = Spectrum_Graph(Spectrum, Amino_Mass_Dict)

	Spectrum = np.array([0] + Spectrum)

	paths = list(nx.all_simple_paths(G,0,Spectrum[-1]))
	name_paths = [''.join([list(G[item[i-1]][item[i]][0].values())[0] for i in range(1,len(item))]) for item in paths]

	def IdealSpectrum(Peptide):
	    prefs = [Peptide[0:i] for i in range(len(Peptide) + 1)]
	    suffs = [Peptide[::-1][0:i][::-1] for i in range(1,len(Peptide))]
	    both = prefs + suffs
	    ans = [sum(list(map(lambda x: Amino_Mass_Dict[x], item))) for item in both]
	    return(sorted(ans))


	for item in name_paths:
		try:
			if np.all(IdealSpectrum(item) == Spectrum):
				return(item)
		except:
			pass
	

	return("No Ideal Peptide Found")



#usage
print(DecodingIdealSpectrum([57,71,154,185,301,332,415,429,486], AminoAcidMass))







'''
Converting a Peptide into a Peptide Vector Problem: Convert a peptide into a peptide vector.

Input: An amino acid string Peptide.
Output: The peptide vector Peptide'.
Code Challenge: Solve the Converting a Peptide into a Peptide Vector Problem.

Input: An amino acid string P.
Output: The peptide vector of P (in the form of space-separated integers).
'''


def Peptide_to_Vector(Peptide, Amino_Mass_Dict):

	prefixes = [Peptide[0:i] for i in range(len(Peptide) + 1)]
	ans = sorted([sum(list(map(lambda x: Amino_Mass_Dict[x], item))) for item in prefixes])
	one_hot = [0] * ans[-1]

	for item in ans:
		one_hot[item-1] = 1

	return(one_hot)


#usage
print(' '.join([str(item) for item in (Peptide_to_Vector('XZZXX', minipeptides))]))


'''
Converting a Peptide Vector into a Peptide Problem: Convert a peptide vector into a peptide.

Input: A binary vector P.
Output: A peptide whose peptide vector is equal to P (if such a peptide exists).
Code Challenge: Solve the Converting a Peptide Vector into a Peptide Problem.

Input: A space-delimited binary vector P.
Output: An amino acid string whose binary peptide vector matches P. For masses with more than one amino acid, any choice may be used.
'''


def Vector_to_Peptide(Peptide_Vector, Amino_Mass_Dict):
	
	revAminoAcidMass = {item[1] : item[0] for item in Amino_Mass_Dict.items()} 

	
	Spectrum = [i+1 for i, item in enumerate(Peptide_Vector) if item == 1]
	
	peptide = revAminoAcidMass[Spectrum[0]]

	for i in range(1,len(Spectrum)):
		peptide += revAminoAcidMass[Spectrum[i] - Spectrum[i-1]]

	return(peptide)


#usage
print(Vector_to_Peptide([0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1], minipeptides))






'''
Code Challenge: Solve the Peptide Sequencing Problem.

Input: A space-delimited spectral vector Spectrum'.
Output: An amino acid string with maximum score against Spectrum'. For masses with more than one amino acid, any choice may be used.
'''


def PeptideSequencing(Spectrum, Amino_Mass_Dict, Manual = False):

	revAminoAcidMass = {item[1] : item[0] for item in Amino_Mass_Dict.items()} 
	Spectrum = [0] + Spectrum
	indices = [i for i in range(len(Spectrum))]
	outer_prod = np.tril(np.subtract.outer(indices,indices))

	if Manual == False:
		G = nx.DiGraph()
		

		for i in range(len(Spectrum)):
		    G.add_node(i, weight = 0)

	
		for i in range(outer_prod.shape[0]):
		    for j in range(outer_prod.shape[0]):
		        if outer_prod[i,j] in revAminoAcidMass.keys():
		            G.add_edge(j,i)

		for i in range(len(Spectrum)):
		    G.nodes[i]['weight'] = -Spectrum[i]


		bestpath = (nx.bellman_ford_path(G,0,len(G.nodes)-1, weight = lambda x,y,z : G.nodes[y]['weight']))

		bestpeptide = [revAminoAcidMass[bestpath[i+1] - bestpath[i]] for i in range(len(bestpath) - 1)]

		return(bestpeptide)

	else:

		l = len(Spectrum)
		adj = [[] for _ in range(l)]
		for i in range(outer_prod.shape[0]):
		    for j in range(outer_prod.shape[0]):
		        if outer_prod[i,j] in revAminoAcidMass.keys():
		            adj[j].append(i)

		# Bellman-Ford 
		dist = [-float('inf')] * l
		parent = [None] * l
		dist[0] = 0
		updated = True
		for i in range(l-1):
		    if not updated:
		        break
		    updated = False
		    for u in range(l):
		        for v in adj[u]:
		            if dist[u] + Spectrum[v] > dist[v]:
		                dist[v] = dist[u] + Spectrum[v]
		                parent[v] = u
		                updated = True


		u = l-1
		bestpath = [u]
		while 0 != u:
		    u = parent[u]
		    bestpath = [u] + bestpath

		bestpeptide = [revAminoAcidMass[bestpath[i+1] - bestpath[i]] for i in range(len(bestpath) - 1)]
		return(bestpeptide)







#usage
testcase = '0 0 0 4 -2 -3 -1 -7 6 5 3 2 1 9 3 -8 0 3 1 2 1 8'.split(' ')
testcase = [int(item) for item in testcase]
print(''.join(PeptideSequencing(testcase, minipeptides, Manual = True)))





