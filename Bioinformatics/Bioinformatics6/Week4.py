import numpy as np 
import pandas as pd 
from functools import reduce
import string

#The algorithms here are slower because of the use of Pandas but I think it increases readability



'''
Helper Functions for loading data
loading data for problems 1 and 2
'''
def load_matrix_data_prob_1_2(filepath):

	with open(filepath) as f:
		lines = (f.readlines())


	lineseps = [i for i, item in enumerate(lines) if '-' in item]
	matrix_data = lines[int(lineseps[-1]) + 1:]
	for i in range(len(matrix_data[1:])):
		matrix_data[i+1] = matrix_data[i+1][:-1]
	
	columns = matrix_data[0].split('\t')[1:]
	columns[-1] = columns[-1][:-1]
	rows = [item[0] for item in matrix_data[1:]]
	
	data = [[float(z) for z in item.split('\t')[1:] if z not in string.printable] for item in matrix_data[1:]]

	return(pd.DataFrame(np.array(data), rows, columns))



'''
Probability of a Hidden Path Problem: Compute the probability of a hidden path.

Input: A hidden path π in an HMM (Σ, States, Transition, Emission).
Output: The probability of this path, Pr(π).
Code Challenge: Solve the Probability of a Hidden Path Problem.

Input: A hidden path π followed by the states States and transition matrix Transition of an HMM (Σ, States, Transition, Emission).
Output: The probability of this path, Pr(π).
Note: You may assume that transitions from the initial state occur with equal probability.

'''


def Prob_Hidden_Path(filepath):
	df = load_matrix_data_prob_1_2(filepath)
	initprob = 1/df.shape[0]
	path = open(filepath).readlines()[0]
	preans = reduce(lambda x, y: x * y, (df.loc[path[i-1], path[i]] for i in range(1,len(path) - 1)))
	return(preans * initprob)







'''
Probability of an Outcome Given a Hidden Path Problem: Compute the probability that an HMM will emit a string given its hidden path.

Input: A string x = x1 . . . xn emitted by an HMM (Σ, States, Transition, Emission) and a hidden path π = π1 . . . πn.
Output: The conditional probability Pr(x|π) that x will be emitted given that the HMM follows the hidden path π.
Code Challenge: Solve the Probability of an Outcome Given a Hidden Path Problem.

Input: A string x, followed by the alphabet from which x was constructed, followed by a hidden path π, followed by the states States and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
Output: The conditional probability Pr(x|π) that x will be emitted given that the HMM follows the hidden path π.
Note: You may assume that transitions from the initial state occur with equal probability.
'''


def Prob_Outcome_Given_Path(filepath):
	df = load_matrix_data_prob_1_2(filepath)
	path = open(filepath).readlines()[4]
	emission = open(filepath).readlines()[0]
	ans = reduce(lambda x, y: x * y, (df.loc[path[i],emission[i]] for i in range(len(path) - 1)))
	return ans





'''
Helper Functions for loading data
loading data for problems 3 and 4
'''

	

def load_matrix_data_prob_3_4(filepath):

	with open(filepath) as f:
		lines = (f.readlines())
	path = lines[0][:-1]
	lineseps = [i for i, item in enumerate(lines) if '-' in item]
	emissions = lines[lineseps[0] + 1].split(' ')
	emissions[-1] = emissions[-1][:-1]
	States = lines[lineseps[1] + 1].split(' ')
	States[-1] = States[-1][:-1]
	mat1 = lines[lineseps[2] + 1:lineseps[3]][1:]
	mat2 = lines[lineseps[3] + 1:][1:]
	mat1 = [[float(mixed) for mixed in item.split('\t') if mixed not in string.printable] for item in mat1]
	mat2 = [[float(mixed) for mixed in item.split('\t') if mixed not in string.printable] for item in mat2]


	TransitionMatrix = pd.DataFrame(mat1, States,States)
	EmissionMatrix = pd.DataFrame(mat2, States, emissions)

	return(path, ''.join(States), TransitionMatrix, EmissionMatrix)





'''
Code Challenge: Implement the Viterbi algorithm solving the Decoding Problem.

Input: A string x, followed by the alphabet from which x was constructed, followed by the states States, transition matrix Transition, and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
Output: A path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π.
Note: You may assume that transitions from the initial state occur with equal probability.
'''


#AND

'''
Code Challenge: Solve the Outcome Likelihood Problem.

Input: A string x, followed by the alphabet from which x was constructed, followed by the states States, transition matrix Transition, and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
Output: The probability Pr(x) that the HMM emits x.
Note: You may assume that transitions from the initial state occur with equal probability.
'''


#Set Likelihood to True if you want the likelihood of the state-path outcome



def Viterbi_Algo(path, states, TransitionMatrix, EmissionMatrix, likelihood = False):
	s_source = 1
	pl = len(path)
	graph_matrix = pd.DataFrame(np.zeros((len(states), len(path))), list(states), list(range(len(path))))
	backtrack = graph_matrix.copy()
	graph_matrix.iloc[:,0] = [(1/len(states))*EmissionMatrix.loc[k][path[0]] for k in states]

	if likelihood == True:
		max_likelihood_matrix = graph_matrix.copy()

	for i in range(1,pl):
		for state in states:
			args = [graph_matrix.loc[l][i-1] * TransitionMatrix.loc[l][state] * EmissionMatrix.loc[state][path[i]] for l in states]
			
			if likelihood == True:
				args2 = [max_likelihood_matrix.loc[l][i-1] * TransitionMatrix.loc[l][state] * EmissionMatrix.loc[state][path[i]] for l in states]


			graph_matrix.loc[state][i] = max(args)
			backtrack.loc[state][i] = np.argmax(args)
			if likelihood == True:
				max_likelihood_matrix.loc[state][i] = sum(args2)


	print(graph_matrix, 'testprob')
	final_state = states[np.argmax(graph_matrix.iloc[:,-1])]
	init = backtrack.loc[final_state][len(path)-1]
	print(init)
	ans = [final_state, states[int(init)]]
	for i in range(len(path) - 2,0,-1):
		state = states[int(init)]
		init = backtrack.loc[state][i]
		ans.append(states[int(init)])

	if likelihood == True:
		print(max_likelihood_matrix)
		return(ans[::-1], sum(max_likelihood_matrix.iloc[:,-1]))
	else:
		return(ans[::-1])


