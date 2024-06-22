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


#----------------------------------------------

'''
Protein Translation Problem: Translate an RNA string into an amino acid string.

Input: An RNA string Pattern and the array GeneticCode.
Output: The translation of Pattern into an amino acid string Peptide.
'''

def protein_translation(rna_pattern):
    protein = ""
    i = 0
    while i < len(rna_pattern):
        codon = rna_pattern[i:i + 3]
        amino_acid = rna_to_protein[codon]
        if amino_acid == "Stop":
            break
        protein += amino_acid
        i += 3
    return protein


'''
Helper Function:
Reverse compliment
'''
def reverse_compliment(inp):
	complimentdict = {'A' : 'U', "C" : "G", "G" : "C", 'U' : 'A'}
	return("".join([complimentdict[x] for x in inp][::-1]))

'''
Peptide Encoding Problem: Find substrings of a genome encoding a given amino acid sequence.

Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
Output: All substrings of Text encoding Peptide (if any such substrings exist).
Code Challenge: Solve the Peptide Encoding Problem. Click here for the RNA codon table corresponding to the array GeneticCode.
'''

def peptide_encoding(dna, peptide):
    rna = dna.replace("T", "U")
    peptide_len = len(peptide)
    substrings = []
    for i in range(len(dna) - peptide_len * 3 + 1):
        if protein_translation(rna[i:i + peptide_len * 3]) == peptide:
            substrings.append(dna[i:i + peptide_len * 3])
        if protein_translation(reverse_compliment(rna[i:i + peptide_len * 3])) == peptide:
            substrings.append(dna[i:i + peptide_len * 3])
    return substrings



'''
Exercise Break: How many subpeptides does a cyclic peptide of length n have?

Input: An integer n.
Output: The number of subpeptides of a cyclic peptide of length n.
'''

def count_subpeptides(n):
    return(n * (n-1))


'''
Generating Theoretical Spectrum Problem: Generate the theoretical spectrum of a cyclic peptide.

Input: An amino acid string Peptide.
Output: Cyclospectrum(Peptide).
Code Challenge: Solve the Generating Theoretical Spectrum Problem.
'''


def gen_theoretical_spectrum(peptide):
	n = len(peptide)
	peptide_weight = sum(AminoAcidMass[item] for item in peptide)
	pos_subpeptides = [(peptide * 2)[i : i + j] for i in range(n) for j in range(1,n)]
	weights = sorted([(sum([AminoAcidMass[peptide] for peptide in item])) for item in pos_subpeptides])
	return([0] + weights + [(peptide_weight)])




'''
Counting Peptides with Given Mass Problem: Compute the number of peptides of given mass.

Input: An integer m.
Output: The number of linear peptides having integer mass m.
'''
def CountingPeptides(n, masstable = AminoAcidMass):
    masstable = list(set(masstable.values()))
    m = len(masstable)
    table = [0] * (n+1)
    table[0] = 1
    for i in range(n+1):
        currSum = 0
        for j in range(m):
            if i - masstable[j] >= 0:
                currSum += table[i-masstable[j]]
        table[i] += currSum    
    return table[n]


'''
Exercise Break: How many subpeptides does a linear peptide of given length n have? (Include the empty peptide and the entire peptide.)

Input: An integer n.
Output: The number of subpeptides of a linear peptide of length n.
'''
def CountingSubpeptides(n):
    return(n * (n+1)//2 + 1)


'''
Code Challenge: Implement LinearSpectrum.

Input: An amino acid string Peptide.
Output: The linear spectrum of Peptide.
'''

def LinearSpectrum(Peptide, AminoAcidMass = AminoAcidMass):
    PrefixMass = [0] * (len(Peptide)+1) 
    for i in range(0,len(Peptide)):
        for item in list(AminoAcidMass.keys()):
            if item == Peptide[i]:
              
                PrefixMass[i] = (PrefixMass[i-1] + AminoAcidMass[item])
    
    PrefixMass = sorted(PrefixMass)
    
    LinearSpectrum = [0]
    
    for i in range(0, len(Peptide)-1):
        for j in range(i+1, len(Peptide)+1):
            LinearSpectrum.append(PrefixMass[j] - PrefixMass[i])

    LinearSpectrum.append(AminoAcidMass[Peptide[-1]])

    return(sorted(LinearSpectrum))



'''
Helper Functions:
Expand - 
For every peptide in a list, it adds a copy to the list each with a candidate peptide concatenated to it
Mass - 
Gets the mass of a Peptide
gen_cyclospectrum -
generates the cyclospectrum of a peptide
consistent - 
Checks if a peptide calculated spectrum is consistent with its theoretical spectrum 
'''

def Expand(candidates):

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

def gen_cyclospectrum(peptide):
	n = len(peptide)
	peptide_weight = sum(AminoAcidMass[item] for item in peptide)
	pos_subpeptides = [(peptide * 2)[i : i + j] for i in range(n) for j in range(1,n)]
	weights = sorted([(sum([AminoAcidMass[peptide] for peptide in item])) for item in pos_subpeptides])
	return([0] + weights + [(peptide_weight)])


def consistent(peptide, Spectrum):
    pepspec = LinearSpectrum(peptide)
    return(all([pepspec.count(item) <= Spectrum.count((item)) for item in pepspec]))


'''
Code Challenge: Implement CyclopeptideSequencing.

'''


def CyclopeptideSequencing(Spectrum):
	CandidatePeptides = []
	FinalPeptides = []
	Flag = True
	while len(CandidatePeptides) != 0 or Flag:
		Flag = False 
		CandidatePeptides = Expand(CandidatePeptides)
		deletions = []
		for peptide in CandidatePeptides:

			if Mass(peptide) == Spectrum[-1]:
				if (gen_cyclospectrum(peptide) == Spectrum):
					toadd = [str(AminoAcidMass[item]) for item in peptide]
					if toadd not in FinalPeptides:
						FinalPeptides.append(toadd)
					deletions.append(peptide)
			elif not consistent(peptide, Spectrum):
				deletions.append(peptide)
		for delete in deletions:
			CandidatePeptides.remove(delete)
	return(FinalPeptides)




