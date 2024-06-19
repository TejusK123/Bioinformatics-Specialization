import numpy as np 
from itertools import accumulate
from collections import Counter



'''
Minimum Skew Problem: Find a position in a genome where the skew diagram attains a minimum.

Input: A DNA string Genome.
Output: All integer(s) i minimizing Skewi (Genome) among all values of i (from 0 to |Genome|).
'''

def GC(x,y):
	if y == 'G':
		try:
			return x + 1
		except:
			return 1
	elif y == 'C':
		try:
			return x - 1
		except:
			return -1
	else:
		try:
			return(x + 0)
		except:
			return 0

def min_max_skew(genome, type_skew = 'min'):
	skewlist = (list(accumulate(genome, GC)))
	skewlist[0] = 0
	skewlist = np.array(skewlist)
	if type_skew == 'min':

		minimum = min(skewlist)

		return((list(np.where(skewlist == minimum))[0]))

	elif type_skew == 'max':
		maximum = max(skewlist)

		return((list(np.where(skewlist == minimum))[0]))

	else:
		return("Invalid Parameter (2)")
	

'''
Hamming Distance Problem: Compute the Hamming distance between two strings.

Input: Two strings of equal length.
Output: The Hamming distance between these strings.
'''

def hamming(seq1, seq2):
	return(sum(np.array(list(seq1)) != np.array(list(seq2))))


'''
Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.

Input: Strings Pattern and Text along with an integer d.
Output: All starting positions where Pattern appears as a substring of Text with at most d mismatches.

--------------------------

Code Challenge: Implement ApproximatePatternCount.

Input: Strings Pattern and Text as well as an integer d.
Output: Countd(Text, Pattern).

'''


def approx_kmer(target, sequence, d):
	ans = []
	length = len(target)
	for i in range(len(sequence)- length + 1):
		if hamming(target, sequence[i : i + length]) <= d:
			ans.append(i)

	ans = [str(x) for x in ans]

	return(ans, len(ans))



'''
Code Challenge: Implement Neighbors to find the d-neighborhood of a string.

Input: A string Pattern and an integer d.
Output: The collection of strings Neighbors(Pattern, d). (You may return the strings in any order, but each line should contain only one string.)
'''

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
			if hamming(pattern[1:], item) < d:
				for nucleotide in nucleotides:
					neighbors.add(nucleotide + item)
			else:
				neighbors.add(pattern[0] + item)

		return(neighbors)

'''
Frequent Words with Mismatches and Reverse Complements Problem: Find the most frequent k-mers (with mismatches and reverse complements) in a string.

Input: A DNA string Text as well as integers k and d.
Output: All k-mers Pattern maximizing the sum Countd(Text, Pattern)+ Countd(Text, Patternrc) over all possible k-mers.
'''

def reverse_compliment(inp):
	complimentdict = {'A' : 'T', "T" : 'A', "C" : "G", "G" : "C"}
	return("".join([complimentdict[x] for x in inp][::-1]))

def frequentkmers_w_mismatches_and_reversecomplement(seq, k, d):
	target_patterns = ""
	cur_map = Counter({})

	for i in range(len(seq) - k):
		kmer = seq[i: i + k]
		neighbors = Neighbors(kmer, d)
		for item in neighbors:
			cur_map[item] += 1

	counts = [(cur_map[item] + cur_map[reverse_compliment(item)], item) for item in cur_map.keys()]
	print(max(counts))
	ans = []
	maxnum = max(counts)[0]
	for item in counts:
		if item[0] == maxnum:
			ans.append(item[1])

	return(' '.join(ans))
