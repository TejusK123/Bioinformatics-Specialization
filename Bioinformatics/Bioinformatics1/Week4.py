import numpy as np 
from scipy.stats import entropy
import random 



'''
Helper Functions
'''

def Motifs(Profile, Dna):
	ans = []
	print(Profile.shape[1])
	for i in range(len(Dna)):
		ans.append(profile_most_probable(Dna[i], Profile.shape[1], Profile))

	return(ans)

def Random_Kmers(Dna, k):
	result = []
	for item in Dna:
		random_i = random.randint(0, len(item) - k)
		result.append(item[random_i : random_i + k])

	return(result)


def Score(motifs, typescore = 'basic'):
	#Basic Score: Number of nucleotides deviating from most common
	#Entropy: Entropy of nucleotides

	if typescore == 'basic':
	    total_score = 0
	    k = len(motifs[0])
	    t = len(motifs)
	    
	    for i in range(k):
	        column = [motif[i] for motif in motifs]
	        max_count = max(column.count(x) for x in set(column))
	        total_score += t - max_count 
	    return total_score


def consensus_profile(Dna):

	motif_matrix = np.stack([np.array(list(item)) for item in Dna])
	profile = np.zeros((4,motif_matrix.shape[1]))

	for i, base in enumerate('ACGT'):
		for j in range(motif_matrix.shape[1]):
			profile[i,j] = (sum(motif_matrix[:,j] == base) + 1)/(motif_matrix.shape[0] + 4)
	return(profile)

'''
Code Challenge: Implement RandomizedMotifSearch.

Input: Integers k and t, followed by a space-separated collection of strings Dna.
Output: A collection BestMotifs resulting from running RandomizedMotifSearch(Dna, k, t) 1,000 times. Remember to use pseudocounts!
'''	

def RandomizedMotifSearch(Dna, k, t):
	motifs = []
	for item in Dna:
		rand_i = random.randint(0,len(item) - k)
		motifs.append(item[rand_i : rand_i + k])

	bestmotifs = motifs 

	while True:
		
		Profile = consensus_profile(motifs)
		#print(Profile.shape)
		#motifs = Motifs(Profile, Dna)
		if Score(motifs, type_score = 'entropy') < Score(bestmotifs, type_score = 'entropy'):
			bestmotifs = motifs
		else:
			return(bestmotifs)


def profile_random_kmer(text, k, profile):
    nucleotide_map = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}
    max_probability = -1
    most_probable_kmer = ""
    debug = []
    probabilities = []
    kmers = []
    for i in range(len(text) - k + 1):

        kmer = text[i:i+k]

        probability = 1.0
        for j, nucleotide in enumerate(kmer):
            probability *= profile[nucleotide_map[nucleotide]][j]
        debug.append((probability, kmer))
        probabilities.append(probability)
        kmers.append(kmer)
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer
    #print((debug))
    probabilities = np.array(probabilities)
    normalize_probabilities = probabilities/np.sum(probabilities)
    return(np.random.choice(kmers, p = normalize_probabilities))
    

'''
Code Challenge: Implement GibbsSampler.

Input: Integers k, t, and N, followed by a space-separated collection of strings Dna.
Output: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts. Remember to use pseudocounts!
'''
def Gibbs_Sampler(Dna, k, t, N, type_score = 'basic'):
	motifs = []
	for item in Dna:
		rand_i = random.randint(0,len(item) - k)
		motifs.append(item[rand_i : rand_i + k])

	BestMotifs = motifs
	for j in range(1,N):
		i = np.random.randint(0,t)
		Profile = consensus_profile(motifs[:i] + motifs[i + 1:])
		motifs[i] = profile_random_kmer(Dna[i], k, Profile)
		if Score(motifs, type_score) < Score(BestMotifs, type_score):
			BestMotifs = motifs 

	return(BestMotifs)