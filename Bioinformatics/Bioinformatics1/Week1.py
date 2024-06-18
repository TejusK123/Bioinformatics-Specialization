import re


'''
Code Challenge: Implement PatternCount (reproduced below).
     Input: Strings Text and Pattern.
     Output: Count(Text, Pattern).
'''

def num_kmer(inp, target):
	return(len(re.findall(target, inp)))


'''
Code Challenge: Solve the Frequent Words Problem.

Input: A string Text and an integer k.
Output: All most frequent k-mers in Text.
'''

def most_freq_kmer(inp, k, t = None):
	cur_dict = {}
	for i in range(len(inp) - k + 1):
		kmer = inp[i : i + k]
		if kmer not in cur_dict:
			cur_dict[kmer] = num_kmer(inp, kmer)
	if t == None:
		return(" ".join(list(filter(lambda key: cur_dict[key] == max(cur_dict.values()), cur_dict.keys()))))
	else:
		return(" ".join(list(filter(lambda key: cur_dict[key] >= t, cur_dict.keys()))))
	
'''
Reverse Complement Problem: Find the reverse complement of a DNA string.

Input: A DNA string Pattern.
Output: Patternrc , the reverse complement of Pattern.
'''

def reverse_compliment(inp):
	complimentdict = {'A' : 'T', "T" : 'A', "C" : "G", "G" : "C"}
	return("".join([complimentdict[x] for x in inp][::-1]))



'''
Code Challenge: Solve the Pattern Matching Problem.

Input: Two strings, Pattern and Genome.
Output: A collection of space-separated integers specifying all starting positions where Pattern appears as a substring of Genome.
'''

def Pattern_Matching(pattern, genome):

	#re.finditer only finds non-overlapping patterns 
	#the ?={0} invocation makes it so that finditer does not consume the characters when matching

	ans =  [str(m.start()) for m in re.finditer('(?={0})'.format(pattern), genome)]

	return(" ".join(ans))



'''
Clump Finding Problem: Find patterns forming clumps in a string.

Input: A string Genome, and integers k, L, and t.
Output: All distinct k-mers forming (L, t)-clumps in Genome.
'''

def clump_finding(genome, k, L, t):
	ans = set()
	for i in range(0, len(genome)):
		ans.add(most_freq_kmer(genome[i : i + L], k, t))
	ans.remove('')
	return(" ".join(ans))

print(clump_finding('CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA', 5, 50, 4))
