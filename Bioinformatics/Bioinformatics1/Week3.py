import numpy as np 
from scipy.stats import entropy

'''
Helper Functions:
hamming_distance: compute hamming distance

generate_all_kmers:
Generate all possible k-mers of size k

Neighbors:
Generate all Neighbors of a DNA string
(strings that differ by at most d)

Approx_kmer:
Approximate kmers
'''
def hamming_distance(seq1, seq2):
	return(sum(np.array(list(seq1)) != np.array(list(seq2))))



def generate_all_kmers(nucleotide_options, k):

	permutations = [nucleotide for nucleotide in nucleotide_options]

	for i in range(k - 1):
		new_permutations = []
		for permutation in permutations:
			for nucleotide in nucleotide_options:
				new_permutations.append(permutation + nucleotide)
		permutations = new_permutations

	return(permutations)



def Neighbors(pattern, d):
	nucleotides = {'A', 'T', 'C', 'G'}
	if d == 0:
		return(pattern)

	if len(pattern) == 1:
		return(nucleotides)
	else:
		neighbors = set()
		suffixneighbors = Neighbors(pattern[1:], d)
		for item in suffixneighbors:
			if hamming_distance(pattern[1:], item) < d:
				for nucleotide in nucleotides:
					neighbors.add(nucleotide + item)
			else:
				neighbors.add(pattern[0] + item)

		return(neighbors)


def approx_kmer(target, sequence, d):
	ans = []
	length = len(target)
	for i in range(len(sequence)- length + 1):
		if hamming_distance(target, sequence[i : i + length]) <= d:
			ans.append(i)

	ans = [str(x) for x in ans]

	return(len(ans))

'''
Code Challenge: Implement MotifEnumeration (reproduced below).

Input: Integers k and d, followed by a space-separated collection of strings Dna.
Output: All (k, d)-motifs in Dna.
'''

def MotifEnumeration(DNA, k, d):
	patterns = set()
	for item in DNA: #for each sequence that should contain the motif
		for i in range(len(item) - k + 1): #for each kmer in the sequence
			for j in Neighbors(item[i : i + k], d): #all kmer neighbors which could potentially be the OG kmers just mutated slightly with hamming differences d
				z = [approx_kmer(j, item, d) for item in DNA] #see if there is an approximate match for the kmers in all of the DNA sequences
				if all(z): #If all of the sequences (potentially could change this to some percentage) contain the motif then that is a potential locus of interest
					patterns.add(j)

	return(' '.join(patterns))


'''
Code Challenge: Implement DistanceBetweenPatternAndStrings.

Input: A string Pattern followed by a collection of space-separated strings Dna.
Output: d(Pattern, Dna).
'''

def DistanceBetweenPatternAndStrings(Pattern, Dna):
	k = len(Pattern)
	sum_distance = 0

	for sequence in Dna:
		cur_distance = float('inf')
		for i in range(len(sequence) - k + 1):
			cur_distance = min(hamming_distance(sequence[i : i + k], Pattern), cur_distance)
		sum_distance += cur_distance

	return(sum_distance)

'''
Code Challenge: Implement MedianString.

Input: An integer k, followed by a space-separated collection of strings Dna.
Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers. (If there are multiple such strings Pattern, then you may return any one.)
'''

def MedianString(Dna, k):

	potential_motifs = []

	max_dist = float('inf')

	for kmer in generate_all_kmers('ACTG', k):

		kmer_distance = DistanceBetweenPatternAndStrings(kmer, Dna)

		if kmer_distance <= max_dist:
			potential_motifs.append((kmer_distance, kmer))
			max_dist = kmer_distance

	return(potential_motifs)






'''
Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.

Input: A string Text, an integer k, and a 4 × k matrix Profile.
Output: A Profile-most probable k-mer in Text.
'''
def profile_most_probable(text, k, profile):
    nucleotide_map = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}
    max_probability = -1
    most_probable_kmer = ""
    debug = []
    for i in range(len(text) - k + 1):

        kmer = text[i:i+k]

        probability = 1.0
        for j, nucleotide in enumerate(kmer):
            probability *= profile[nucleotide_map[nucleotide]][j]
        debug.append((probability, kmer, profile))
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer

    return most_probable_kmer



'''
Helper Function:
consensus_profile: Generate a consensus profile given a DNA string

Score: Generate the Score of a motif
'''

def consensus_profile(Dna):

	motif_matrix = np.stack([np.array(list(item)) for item in Dna])
	profile = np.zeros((4,motif_matrix.shape[1]))

	for i, base in enumerate('ACGT'):
		for j in range(motif_matrix.shape[1]):
			#Delete + 4 for no pseudocounts
			profile[i,j] = (sum(motif_matrix[:,j] == base) + 1)/(motif_matrix.shape[0] + 4)
	return(profile)

def Score(motifs, type_score = 'entropy'):
	if type_score == 'entropy':
		return(sum(entropy(consensus_profile(motifs))))
	if type_score == 'basic':
		total_score = 0
		k = len(motifs[0])
		t = len(motifs)

		for i in range(k):
			column = [motif[i] for motif in motifs]
			max_count = max(column.count(x) for x in set(column))
			total_score += t - max_count 
		return total_score

'''
Code Challenge: Implement GreedyMotifSearch with pseudocounts.

Input: Integers k and t, followed by a space-separated collection of strings Dna.
Output: A collection of strings BestMotifs resulting from applying GreedyMotifSearch(Dna, k, t) with pseudocounts. If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.
'''

def greedy_motif_search_with_pseudocounts(Dna, k, t):
	BestMotifs = [seq[:k] for seq in Dna]
	for i in range(len(Dna[0]) - k + 1):
		Motifs = [Dna[0][i : i + k]]
		for j in range(1,t):
			Profile = consensus_profile(Motifs)
			Motifs.append(profile_most_probable(Dna[j], k, Profile))

		if Score(Motifs) < Score(BestMotifs):
			BestMotifs = Motifs 

	return(BestMotifs)















