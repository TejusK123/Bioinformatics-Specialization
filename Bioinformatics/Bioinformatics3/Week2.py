import numpy as np

#Scoring Matrices:
PAM250 = {'A': {'A': 2, 'C': -2, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -2, 'N': 0, 'Q': 0, 'P': 1, 'S': 1, 'R': -2, 'T': 1, 'W': -6, 'V': 0, 'Y': -3},
'C': {'A': -2, 'C': 12, 'E': -5, 'D': -5, 'G': -3, 'F': -4, 'I': -2, 'H': -3, 'K': -5, 'M': -5, 'L': -6, 'N': -4, 'Q': -5, 'P': -3, 'S': 0, 'R': -4, 'T': -2, 'W': -8, 'V': -2, 'Y': 0}, 
'E': {'A': 0, 'C': -5, 'E': 4, 'D': 3, 'G': 0, 'F': -5, 'I': -2, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 
'D': {'A': 0, 'C': -5, 'E': 3, 'D': 4, 'G': 1, 'F': -6, 'I': -2, 'H': 1, 'K': 0, 'M': -3, 'L': -4, 'N': 2, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 
'G': {'A': 1, 'C': -3, 'E': 0, 'D': 1, 'G': 5, 'F': -5, 'I': -3, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -3, 'T': 0, 'W': -7, 'V': -1, 'Y': -5}, 
'F': {'A': -3, 'C': -4, 'E': -5, 'D': -6, 'G': -5, 'F': 9, 'I': 1, 'H': -2, 'K': -5, 'M': 0, 'L': 2, 'N': -3, 'Q': -5, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -1, 'Y': 7}, 
'I': {'A': -1, 'C': -2, 'E': -2, 'D': -2, 'G': -3, 'F': 1, 'I': 5, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -2, 'S': -1, 'R': -2, 'T': 0, 'W': -5, 'V': 4, 'Y': -1}, 
'H': {'A': -1, 'C': -3, 'E': 1, 'D': 1, 'G': -2, 'F': -2, 'I': -2, 'H': 6, 'K': 0, 'M': -2, 'L': -2, 'N': 2, 'Q': 3, 'P': 0, 'S': -1, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': 0}, 
'K': {'A': -1, 'C': -5, 'E': 0, 'D': 0, 'G': -2, 'F': -5, 'I': -2, 'H': 0, 'K': 5, 'M': 0, 'L': -3, 'N': 1, 'Q': 1, 'P': -1, 'S': 0, 'R': 3, 'T': 0, 'W': -3, 'V': -2, 'Y': -4}, 
'M': {'A': -1, 'C': -5, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 2, 'H': -2, 'K': 0, 'M': 6, 'L': 4, 'N': -2, 'Q': -1, 'P': -2, 'S': -2, 'R': 0, 'T': -1, 'W': -4, 'V': 2, 'Y': -2}, 
'L': {'A': -2, 'C': -6, 'E': -3, 'D': -4, 'G': -4, 'F': 2, 'I': 2, 'H': -2, 'K': -3, 'M': 4, 'L': 6, 'N': -3, 'Q': -2, 'P': -3, 'S': -3, 'R': -3, 'T': -2, 'W': -2, 'V': 2, 'Y': -1}, 
'N': {'A': 0, 'C': -4, 'E': 1, 'D': 2, 'G': 0, 'F': -3, 'I': -2, 'H': 2, 'K': 1, 'M': -2, 'L': -3, 'N': 2, 'Q': 1, 'P': 0, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -2, 'Y': -2}, 
'Q': {'A': 0, 'C': -5, 'E': 2, 'D': 2, 'G': -1, 'F': -5, 'I': -2, 'H': 3, 'K': 1, 'M': -1, 'L': -2, 'N': 1, 'Q': 4, 'P': 0, 'S': -1, 'R': 1, 'T': -1, 'W': -5, 'V': -2, 'Y': -4}, 
'P': {'A': 1, 'C': -3, 'E': -1, 'D': -1, 'G': 0, 'F': -5, 'I': -2, 'H': 0, 'K': -1, 'M': -2, 'L': -3, 'N': 0, 'Q': 0, 'P': 6, 'S': 1, 'R': 0, 'T': 0, 'W': -6, 'V': -1, 'Y': -5}, 
'S': {'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': -1, 'P': 1, 'S': 2, 'R': 0, 'T': 1, 'W': -2, 'V': -1, 'Y': -3}, 
'R': {'A': -2, 'C': -4, 'E': -1, 'D': -1, 'G': -3, 'F': -4, 'I': -2, 'H': 2, 'K': 3, 'M': 0, 'L': -3, 'N': 0, 'Q': 1, 'P': 0, 'S': 0, 'R': 6, 'T': -1, 'W': 2, 'V': -2, 'Y': -4}, 
'T': {'A': 1, 'C': -2, 'E': 0, 'D': 0, 'G': 0, 'F': -3, 'I': 0, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -1, 'T': 3, 'W': -5, 'V': 0, 'Y': -3}, 
'W': {'A': -6, 'C': -8, 'E': -7, 'D': -7, 'G': -7, 'F': 0, 'I': -5, 'H': -3, 'K': -3, 'M': -4, 'L': -2, 'N': -4, 'Q': -5, 'P': -6, 'S': -2, 'R': 2, 'T': -5, 'W': 17, 'V': -6, 'Y': 0}, 
'V': {'A': 0, 'C': -2, 'E': -2, 'D': -2, 'G': -1, 'F': -1, 'I': 4, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -1, 'S': -1, 'R': -2, 'T': 0, 'W': -6, 'V': 4, 'Y': -2}, 
'Y': {'A': -3, 'C': 0, 'E': -4, 'D': -4, 'G': -5, 'F': 7, 'I': -1, 'H': 0, 'K': -4, 'M': -2, 'L': -1, 'N': -2, 'Q': -4, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -2, 'Y': 10}
}


BLOSUM62 = {'A': {'A': 4, 'C': 0, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0, 'W': -3, 'V': 0, 'Y': -2}, 
'C': {'A': 0, 'C': 9, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3, 'K': -3, 'M': -1, 'L': -1, 'N': -3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2}, 
'E': {'A': -1, 'C': -4, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2}, 
'D': {'A': -2, 'C': -3, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3, 'H': -1, 'K': -1, 'M': -3, 'L': -4, 'N': 1, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -4, 'V': -3, 'Y': -3}, 
'G': {'A': 0, 'C': -3, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3}, 
'F': {'A': -2, 'C': -2, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0, 'H': -1, 'K': -3, 'M': 0, 'L': 0, 'N': -3, 'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3}, 
'I': {'A': -1, 'C': -1, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4, 'H': -3, 'K': -3, 'M': 1, 'L': 2, 'N': -3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1}, 
'H': {'A': -2, 'C': -3, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H': 8, 'K': -1, 'M': -2, 'L': -3, 'N': 1, 'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2}, 
'K': {'A': -1, 'C': -3, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1, 'K': 5, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': -2}, 
'M': {'A': -1, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1, 'H': -2, 'K': -1, 'M': 5, 'L': 2, 'N': -2, 'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': 1, 'Y': -1}, 
'L': {'A': -1, 'C': -1, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2, 'H': -3, 'K': -2, 'M': 2, 'L': 4, 'N': -3, 'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -2, 'V': 1, 'Y': -1}, 
'N': {'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 6, 'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2}, 
'Q': {'A': -1, 'C': -3, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0, 'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2, 'V': -2, 'Y': -1}, 
'P': {'A': -1, 'C': -3, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -2, 'L': -3, 'N': -2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T': -1, 'W': -4, 'V': -2, 'Y': -3}, 
'S': {'A': 1, 'C': -1, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1, 'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2}, 
'R': {'A': -1, 'C': -3, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 2, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2}, 
'T': {'A': 0, 'C': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5, 'W': -2, 'V': 0, 'Y': -2}, 
'W': {'A': -3, 'C': -2, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3, 'H': -2, 'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2}, 
'V': {'A': 0, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3, 'H': -3, 'K': -2, 'M': 1, 'L': 1, 'N': -3, 'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1}, 
'Y': {'A': -2, 'C': -2, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1, 'H': 2, 'K': -2, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7}
}





'''
Code Challenge: Solve the Global Alignment Problem.

Input: A match reward, a mismatch penalty, an indel penalty, and two nucleotide strings.
Output: The maximum alignment score of these strings followed by an alignment achieving this maximum score.
'''

def Global_Alignment(v, w, match_reward, mu, sigma):

	def LCSBackTrack(v, w, match_reward, mu, sigma):
		n = len(v)
		m = len(w)
		s = np.zeros((n+1,m+1))
		s[:,0] = np.array(list(range(0, -sigma * (n+1), -sigma)))
		s[0,:] = np.array(list(range(0, -sigma * (m+1), -sigma)))
		backtrack = np.zeros((n+1,m+1))

		backtrack[:,0] = 1
		backtrack[0,:] = 2
		backtrack[0,0] = 0
		for i in range(1, n+1):
			for j in range(1, m+1):
				match = 0 
				if v[i-1] == w[j-1]:
					match = match_reward
				else:
					match = -mu
				s[i,j] = max(s[i-1,j] - sigma, s[i,j-1] - sigma, s[i-1,j-1] + match)
				if s[i, j] == s[i-1, j-1] + match:
					backtrack[i, j] = 3 #Diag
				elif s[i, j] == s[i-1, j] - sigma:
					backtrack[i, j] = 1 #down
				else:
					backtrack[i,j] = 2


	
		return backtrack,s

	backtrack, s = LCSBackTrack(v,w, match_reward, mu, sigma)
	i = len(v)
	j = len(w)
	s1 = ''
	s2 = ''
	while i > 0 and j > 0:
		if 1 == backtrack[i, j]:
			i -= 1
			s1 += v[i]
			s2 += "-"
		elif 2 == backtrack[i, j]:
			j -= 1
			s1 += "-"
			s2 += w[j]
		else:
			
			i -= 1
			j -= 1
			s1 += v[i]
			s2 += w[j]
			

	return(s1[::-1], s2[::-1], s[-1,-1])



#usage 
t1 = 'GAGA'
t2 = 'GAT'
match_reward = 1
mismatch_penalty = 1
indel_penalty = 2

print(Global_Alignment(t1,t2,match_reward, mismatch_penalty, indel_penalty))






'''
Code Challenge: Solve the Local Alignment Problem.

Input: Two protein strings written in the single-letter amino acid alphabet.
Output: The maximum score of a local alignment of the strings, followed by a local alignment of these strings achieving the maximum score. Use the PAM250 scoring matrix for matches and mismatches as well as the indel penalty σ = 5.
'''




def Local_Alignment(v, w, scoring_matrix, sigma):


	def LCSBackTrack(v, w, scoring_matrix, sigma):
		n = len(v)
		m = len(w)
		s = np.zeros((n+1,m+1))
		backtrack = np.zeros((n+1,m+1))
		for i in range(1, n+1):
			for j in range(1, m+1):
				match = 0 
				if v[i-1] == w[j-1]:
					match = scoring_matrix[v[i-1]][w[j-1]]
				else:
					match = scoring_matrix[v[i-1]][w[j-1]]
				s[i,j] = max(s[i-1,j] - sigma, s[i,j-1] - sigma, s[i-1,j-1] + match, 0)
				if s[i,j] == 0:
					backtrack[i,j] = 0
				elif s[i,j] == s[i-1,j] - sigma:
					backtrack[i,j] = 1
				elif s[i,j] == s[i,j-1] - sigma:
					backtrack[i,j] = 2
				elif s[i,j] == s[i-1,j-1] + match:
					backtrack[i,j] = 3

		s[n,m] = np.max(s)
		
		return backtrack,s


	backtrack, s = LCSBackTrack(v,w, scoring_matrix, sigma)
	
	i = len(v)
	j = len(w)
	s1 = ''
	s2 = ''
	while i > 0 and j > 0:
		if 1 == backtrack[i, j]:
			i -= 1
			s1 += v[i]
			s2 += "-"
		elif 2 == backtrack[i, j]:
			j -= 1
			s1 += "-"
			s2 += w[j]
		elif 0 == backtrack[i,j]:
			i = 0
			j = 0
		else:
			
			i -= 1
			j -= 1
			s1 += v[i]
			s2 += w[j]
			
			
		#print(s1,s2)

	return(s1[::-1], s2[::-1], s[-1,-1])



#usage 
print(Local_Alignment('MEANLY', 'PENALTY', PAM250, 5))



'''
Edit Distance Problem: Find the edit distance between two strings.

Input: Two strings.
Output: The edit distance between these strings.
Code Challenge: Solve the Edit Distance Problem.
'''

def EditDistance(v, w):
	n = len(v)
	m = len(w)
	s = np.zeros((n+1,m+1))
	backtrack = np.zeros((n+1,m+1))
	s[:,0] = np.array(list(range(n+1)))
	s[0,:] = np.array(list(range(m+1)))
	for i in range(1, n+1):
		for j in range(1, m+1):
			match = 0 
			if v[i-1] != w[j-1]:
				match = 1 
			s[i,j] = min(s[i-1,j] + 1, s[i,j-1] + 1, s[i-1,j-1] + match)
			if s[i, j] == s[i-1, j]:
				backtrack[i, j] = 1 #down
			elif s[i, j] == s[i, j-1]:
				backtrack[i, j] = 2 #right
			elif s[i, j] == s[i-1, j-1] + 1 and v[i-1] == w[j-1]:
				backtrack[i, j] = 3 #Diag
	
	return int(s[-1,-1])


#usage
print(EditDistance('GAGA', 'GAT'))





'''
Code Challenge: Solve the Fitting Alignment Problem.

Input: Two amino acid strings.
Output: A highest-scoring fitting alignment between v and w. Use the BLOSUM62 scoring table and an indel penalty equal to 1.
'''


def Fitting_Alignment(v,w, scoring_matrix, sigma):

	def LCSBackTrack(v, w, scoring_matrix, sigma):
		n = len(v)
		m = len(w)
		s = np.zeros((n+1,m+1))
		backtrack = np.zeros((n+1,m+1))

		s[:, 0] = 0
		s[0,:] = np.array(list(range(0, -sigma * (m+1), -sigma)))


		for i in range(1, n+1):
			for j in range(1, m+1):
				match = 0 
				if v[i-1] == w[j-1]:
					match = scoring_matrix[v[i-1]][w[j-1]]
				else:
					match = scoring_matrix[v[i-1]][w[j-1]]
				s[i,j] = max(s[i-1,j] - sigma, s[i,j-1] - sigma, s[i-1,j-1] + match, 0)
				if s[i,j] == 0:
					backtrack[i,j] = 0
				elif s[i,j] == s[i-1,j] - sigma:
					backtrack[i,j] = 1
				elif s[i,j] == s[i,j-1] - sigma:
					backtrack[i,j] = 2
				elif s[i,j] == s[i-1,j-1] + match:
					backtrack[i,j] = 3

		s[n,m] = np.max(s)
		
		return backtrack,s


	backtrack, s = LCSBackTrack(v,w, scoring_matrix, sigma)
	
	i = np.argmax(s[:,len(w)])
	j = len(w)
	score = s[i,j]

	s1 = ''
	s2 = ''
	while i > 0 and j > 0:
		if 1 == backtrack[i, j]:
			i -= 1
			s1 += v[i]
			s2 += "-"
		elif 2 == backtrack[i, j]:
			j -= 1
			s1 += "-"
			s2 += w[j]
		elif 0 == backtrack[i,j]:
			i = 0
			j = 0
		else:
			
			i -= 1
			j -= 1
			s1 += v[i]
			s2 += w[j]
			
			

	return(s1[::-1], s2[::-1], score)



#usage
print(Fitting_Alignment('DISCREPANTLY', 'PATENT', BLOSUM62, 1))






'''
Code Challenge: Solve the Overlap Alignment Problem.

Input: A match reward, a mismatch penalty, an indel penalty, and two nucleotide strings v and w.
Output: The score of an optimal overlap alignment of v and w, followed by an alignment of a suffix v' of v and a prefix w' of w achieving this maximum score.
'''


def Overlap_Alignment(v, w, match_reward, mu, sigma):

	def LCSBackTrack(v, w, match_reward, mu, sigma):
		n = len(v)
		m = len(w)
		s = np.zeros((n+1,m+1))
		s[0,:] = np.array(list(range(0, -sigma * (m+1), -sigma)))

		backtrack = np.zeros((n+1,m+1))
		backtrack[:,0] = 1
		backtrack[0,:] = 2
		backtrack[0,0] = 0

		for i in range(1, n+1):
			for j in range(1, m+1):
				match = 0 
				if v[i-1] == w[j-1]:
					match = match_reward
				else:
					match = -mu
				s[i,j] = max(s[i-1,j] - sigma, s[i,j-1] - sigma, s[i-1,j-1] + match)
				if s[i, j] == s[i-1, j-1] + match:
					backtrack[i, j] = 3 #Diag
				elif s[i, j] == s[i-1, j] - sigma:
					backtrack[i, j] = 1 #down
				else:
					backtrack[i,j] = 2


		
		return backtrack,s

	backtrack, s = LCSBackTrack(v,w, match_reward, mu, sigma)

	#max not guarenteed to be bottom right corner
	i, j = np.unravel_index(np.argmax(s), s.shape)
	ans = s[i,j]
	s1 = ''
	s2 = ''
	while i > 0 and j > 0:
		if 1 == backtrack[i, j]:
			i -= 1
			s1 += v[i]
			s2 += "-"
		elif 2 == backtrack[i, j]:
			j -= 1
			s1 += "-"
			s2 += w[j]
		else:
			
			i -= 1
			j -= 1
			s1 += v[i]
			s2 += w[j]
			
			
		

	return(s1[::-1], s2[::-1], ans)


#usage
print(Overlap_Alignment('GAGA', 'GAT', 1,1,2))