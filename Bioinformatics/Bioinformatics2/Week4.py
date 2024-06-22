import numpy as np
from collections import Counter
'''
Useful Dictionaries
'''
AminoAcidMass = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113,
                 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
                 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

rna_to_protein = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}


protein_to_rna = {'K': ['AAA', 'AAG'],
'N': ['AAC', 'AAU'],
'T': ['ACA', 'ACC', 'ACG', 'ACU'],
'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGU'],
'S': ['AGC', 'AGU', 'UCA', 'UCC', 'UCG', 'UCU'],
'I': ['AUA', 'AUC', 'AUU'],
'M': ['AUG'],
'Q': ['CAA', 'CAG'],
'H': ['CAC', 'CAU'],
'P': ['CCA', 'CCC', 'CCG', 'CCU'],
'L': ['CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'],
'E': ['GAA', 'GAG'],
'D': ['GAC', 'GAU'],
'A': ['GCA', 'GCC', 'GCG', 'GCU'],
'G': ['GGA', 'GGC', 'GGG', 'GGU'],
'V': ['GUA', 'GUC', 'GUG', 'GUU'],
'Stop': ['UAA', 'UAG', 'UGA'],
'Y': ['UAC', 'UAU'],
'C': ['UGC', 'UGU'],
'W': ['UGG'],
'F': ['UUC', 'UUU']
}

'''
Helper Functions:
gen_cyclospectrum-
generates the cyclospectrum for a protein
LinearSpectrum - 
generates the linearspectrum for a protein
'''

def gen_cyclospectrum(peptide, AminoAcidMass = AminoAcidMass):
	n = len(peptide)
	peptide_weight = sum(AminoAcidMass[item] for item in peptide)
	pos_subpeptides = [(peptide * 2)[i : i + j] for i in range(n) for j in range(1,n)]
	weights = sorted([(sum([AminoAcidMass[peptide] for peptide in item])) for item in pos_subpeptides])
	return([0] + weights + [(peptide_weight)])

def LinearSpectrum(Peptide, AminoAcidMass = AminoAcidMass):
    PrefixMass = [0] * (len(Peptide)+1) 
    for i in range(0,len(Peptide)):
        for item in list(AminoAcidMass.keys()):
            if item == Peptide[i]:
                #print(item)
                PrefixMass[i] = (PrefixMass[i-1] + AminoAcidMass[item])
    
    PrefixMass = sorted(PrefixMass)
    #print(PrefixMass)
    LinearSpectrum = [0]
    
    for i in range(0, len(Peptide)-1):
        for j in range(i+1, len(Peptide)+1):
            LinearSpectrum.append(PrefixMass[j] - PrefixMass[i])

    LinearSpectrum.append(AminoAcidMass[Peptide[-1]])

    return(sorted(LinearSpectrum))

'''
Cyclopeptide Scoring Problem: Compute the score of a cyclic peptide against a spectrum.

Input: An amino acid string Peptide and a collection of integers Spectrum.
Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum).
Code Challenge: Solve the Cyclopeptide Scoring Problem.
'''

def ScoreCyclopeptideSpectra(Peptide, Experimental_Spectra, AminoAcidMass = AminoAcidMass):
	theorectical_spectra = gen_cyclospectrum(Peptide, AminoAcidMass)

	score = 0
	if Peptide == []:
		return(score)
	ts_set = set(theorectical_spectra)
	es_set = set(theorectical_spectra)
	for item in ts_set:
		score += min(theorectical_spectra.count(item), Experimental_Spectra.count(item))

	for item in es_set:
		if item not in ts_set:
			score += Experimental_Spectra.count(item)
	return(score)

'''
Code Challenge: Compute the score of a linear peptide with respect to a spectrum.

Input: An amino acid string Peptide and a collection of integers Spectrum.
Output: The linear score of Peptide with respect to Spectrum, LinearScore(Peptide, Spectrum).
'''

def ScoreLinearCyclopeptideSpectra(Peptide, Experimental_Spectra, AminoAcidMass = AminoAcidMass):
	theorectical_spectra = LinearSpectrum(Peptide, AminoAcidMass)
	if Peptide == []:
		return(0)
	score = 0

	ts_set = set(theorectical_spectra)
	es_set = set(theorectical_spectra)
	for item in ts_set:
		score += min(theorectical_spectra.count(item), Experimental_Spectra.count(item))

	for item in es_set:
		if item not in ts_set:
			score += Experimental_Spectra.count(item)
	return(score)	


'''
Code Challenge: Implement Trim.
Input: A collection of peptides Leaderboard, a collection of integers Spectrum, and an integer N.
Output: The N highest-scoring linear peptides on Leaderboard with respect to Spectrum.
'''


def Trim(Leaderboard, Spectrum, N, AminoAcidMass = AminoAcidMass):

	both = ([(ScoreLinearCyclopeptideSpectra(item, Spectrum, AminoAcidMass), item) for item in Leaderboard])
	
	both = sorted(both, reverse = True)
	leaderboard = [item[1] for item in both]
	scores = [item[0] for item in both]
	

	for j in range(N, len(leaderboard)):
		if scores[j] < scores[N-1]:
			
			leaderboard = leaderboard[:j]

	return(leaderboard)

'''
Helper Functions:
Expand - 
For every peptide in a list, it adds a copy to the list each with a candidate peptide concatenated to it
Mass - 
Gets the mass of a Peptide
extended_mass_table-
generates an extended mass table
'''

def Expand(candidates, AminoAcidMass = AminoAcidMass):

	aminoacids = list(AminoAcidMass.keys())
	if candidates == []:
		return([[item] for item in aminoacids])
	new = []
	for item in candidates:
		for i in range(len(aminoacids)): 
			new.append(item.copy() + [aminoacids[i]])
	return(new)

def Mass(peptide, AminoAcidMass = AminoAcidMass):
	return(sum([AminoAcidMass[item] for item in peptide]))

def extended_mass_table():
    extended_list = {}
    for i in range(57,201):
        extended_list[chr(i)] = str(i)
    return(extended_list)



def LeaderboardCyclopeptideSequencing(Spectrum, N, eight_three = False, AminoAcidMass = AminoAcidMass):
	if eight_three == True:
		eep = []
	Leaderboard = []
	LeaderPeptide = []
	Flag = True
	while len(Leaderboard) != 0 or Flag:
		#print(LeaderPeptide)
		Flag = False 
		Leaderboard = Expand(Leaderboard)
		deletions = []
		for peptide in Leaderboard:
			if eight_three == True:
				if ScoreCyclopeptideSpectra(peptide, Spectrum) == 83:

					eep.append('-'.join([str(AminoAcidMass[item]) for item in peptide]))

			if Mass(peptide) == Spectrum[-1]:
				if ScoreCyclopeptideSpectra(peptide, Spectrum) > ScoreCyclopeptideSpectra(LeaderPeptide, Spectrum):
					LeaderPeptide = peptide
			elif Mass(peptide) > Spectrum[-1]:
				deletions.append(peptide)
		for delete in deletions:
			Leaderboard.remove(delete)

		Leaderboard = Trim(Leaderboard, Spectrum, N)
	LeaderPeptide = [AminoAcidMass[item] for item in LeaderPeptide]
	if eight_three == True:
		print(' '.join(list(set(eep))))
		print(len(list(set(eep))))
	return(LeaderPeptide)





'''
Spectral Convolution Problem: Compute the convolution of a spectrum.

Input: A collection of integers Spectrum in increasing order..
Output: The list of elements in the convolution of Spectrum. If an element has multiplicity k, it should appear exactly k times; you may return the elements in any order.
Code Challenge: Solve the Spectral Convolution Problem.
'''
def spectral_convolution(spectrum):
	spectrum = np.array(spectrum)
	ans = (np.tril(np.subtract.outer(spectrum,spectrum))).flatten()
	return(list(ans[ans != 0]))


'''
Code Challenge: Implement ConvolutionCyclopeptideSequencing.

Input: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.
Output: A cyclic peptide LeaderPeptide with amino acids taken only from the top M elements (and ties) of the convolution of Spectrum that fall between 57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).
'''

def ConvolutionCyclopeptideSequencing(Spectrum, N, M):
	convolutions = Counter(spectral_convolution(Spectrum))
	#print(convolutions)
	z = list(convolutions.items())
	z = (sorted(z, key = lambda x: x[1], reverse = True))
	z = list(filter(lambda x: (x[0] >= 57 and x[0] <= 200), z))
	#print(z)
	aminoacidmasses = {}
	for i in range(M):
		aminoacidmasses[f"aa{z[i][0]}"] = z[i][0]

	last = z[M-1][1]
	i = M-1
	while z[i][1] == last:
		aminoacidmasses[f"aa{z[i][0]}"] = z[i][0]
		i += 1
	
	#print(aminoacidmasses, "\n")
	##############

	Leaderboard = []
	LeaderPeptide = []
	Flag = True
	while len(Leaderboard) != 0 or Flag:
		#print(LeaderPeptide)
		Flag = False 
		Leaderboard = Expand(Leaderboard, aminoacidmasses)
		deletions = []
		for peptide in Leaderboard:
			print(peptide)
			if Mass(peptide, aminoacidmasses) == Spectrum[-1]:
				if ScoreCyclopeptideSpectra(peptide, Spectrum, aminoacidmasses) > ScoreCyclopeptideSpectra(LeaderPeptide, Spectrum, aminoacidmasses):
					LeaderPeptide = peptide
			elif Mass(peptide, aminoacidmasses) > Spectrum[-1]:
				deletions.append(peptide)
		for delete in deletions:
			Leaderboard.remove(delete)

		Leaderboard = Trim(Leaderboard, Spectrum, N, aminoacidmasses)
	LeaderPeptide = [str(aminoacidmasses[item]) for item in LeaderPeptide]
	print('-'.join(LeaderPeptide))
	return(LeaderPeptide)


