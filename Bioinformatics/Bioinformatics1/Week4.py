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


def RandomMotif(Dna, k):
	result = []
	for item in Dna:
		random_i = random.randint(0, len(item) - k)
		result.append(item[random_i : random_i + k])

	return(result)


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

def RandomizedMotifSearch(Dna, k, t, type_score = 'basic'):
	motifs = []
	for item in Dna:
		rand_i = random.randint(0,len(item) - k)
		motifs.append(item[rand_i : rand_i + k])

	bestmotifs = motifs 

	while True:
		
		Profile = consensus_profile(motifs)
		#print(Profile.shape)
		#motifs = Motifs(Profile, Dna)
		if Score(motifs, type_score) < Score(bestmotifs, type_score):
			bestmotifs = motifs
		else:
			return(bestmotifs)



def RepeatedRandomizedMotifSearch(Dna, k, t, type_score = 'basic'):
    best_motifs_overall = " "
    best_score_overall = float('inf')
    for i in range(1000):
        current_motifs = RandomizedMotifSearch(Dna, k, t)
        current_score = Score(current_motifs, type_score = 'basic')
        if current_score < best_score_overall:
            best_motifs_overall = current_motifs
            best_score_overall = current_score
    return best_motifs_overall


'''
Helper Functions
'''
	
def ProfileMostProbableKmer(text, k, profile):
    max_prob = -1
    best_kmer = text[0:k]
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        prob = 1
        for j, nucleotide in enumerate(kmer):
            prob *= profile[nucleotide][j]
            
        if prob > max_prob:
            max_prob = prob
            best_kmer = kmer
    return best_kmer	
	
def Profile(motifs):
    profile_pseudo = {}
    t = len(motifs)
    k = len(motifs[0])
    for symbol in "ACGT":
        profile_pseudo[symbol] = []
        for i in range(k):
            column = [motif[i] for motif in motifs]
            profile_pseudo[symbol].append((column.count(symbol)+1)/ (t+4))
    return profile_pseudo  

'''
Code Challenge: Implement GibbsSampler.

Input: Integers k, t, and N, followed by a space-separated collection of strings Dna.
Output: The strings BestMotifs resulting from running GibbsSampler(Dna, k, t, N) with 20 random starts. Remember to use pseudocounts!
'''
def GibbsSampler(Dna, k, t, N, type_score = 'basic'):
    MotifsRandom = RandomMotif(Dna, k, t)
    BestMotifs = MotifsRandom

    for j in range(N):
        i = random.randint(0, len(Dna) - 1)
        Motifs_no_i = MotifsRandom[:i] + MotifsRandom[i+1:]
        profile = Profile(Motifs_no_i)
        Motif_i = ProfileMostProbableKmer(Dna[i], len(profile['A']), profile)
        MotifsRandom[i] = Motif_i
        if Score(MotifsRandom, type_score) < Score(BestMotifs, type_score):
            BestMotifs = Motifs
    return(BestMotifs)


	
def RepeatedGibbsMotifSearch(Dna, k, t, N, type_score = 'basic'):
    best_motifs_overall = " "
    best_score_overall = float('inf')
    best_entropy_score = float('inf')
    memory = []
    for i in range(20):
        print(i)
        current_motifs = GibbsSampler(Dna, k, t, N, type_score = 'basic')
        current_score = Score(current_motifs)
        entropy_score = sum(entropy(np.array([item for item in Profile(current_motifs).values()])))
        if current_score < best_score_overall:
            best_motifs_overall = current_motifs
            best_score_overall = current_score
            memory.append((current_motifs, current_score, entropy_score))
    print((best_motifs_overall, best_score_overall, entropy_score))
    return(memory)
