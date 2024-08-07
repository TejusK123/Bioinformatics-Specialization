
'''
Helper Function
Outputs Nicely
'''
def output_nicely(Permutation_Sequence):
	eep = '\n'.join([' '.join([str(i) if i < 0 else f'+{i}' for i in item]) for item in Permutation_Sequence])
	return(eep)




'''
Code Challenge: Implement GreedySorting.

Input: A permutation P.
Output: The sequence of permutations corresponding to applying GreedySorting to P, ending with the identity permutation.
'''



def GreedySorting(Permutation, reversal_distance = False):
	approxReversalDistance = 0
	ans = []
	for k in range(len(Permutation)):
		
		if Permutation[k] != (k+1):
			try:
				index = Permutation.index(k+1) 
			except ValueError:
				index = Permutation.index(-(k+1))
			Permutation = Permutation[:k] + [-item for item in Permutation[k:index+1][::-1]] + Permutation[index + 1:]
			
			ans.append(Permutation[:])
			approxReversalDistance += 1
		if Permutation[k] == -(k+1):

			Permutation[k] *= -1
			
			ans.append(Permutation)
			approxReversalDistance += 1
	
	if reversal_distance == True:
		return(ans, approxReversalDistance)
	else:
		return(ans)






'''
Number of Breakpoints Problem: Find the number of breakpoints in a permutation.

Input: A permutation.
Output: The number of breakpoints in this permutation.
Code Challenge: Solve the Number of Breakpoints Problem.
'''

def num_breakpoints(Permutation):
	p1 = 0 
	p2 = 1 
	Permutation = [0] + Permutation + [len(Permutation) + 1]
	num_breakpoints = 0
	while p2 != len(Permutation):
		if Permutation[p2] - Permutation[p1] != 1:
			num_breakpoints += 1
		p1 += 1
		p2 += 1

	return(num_breakpoints)


