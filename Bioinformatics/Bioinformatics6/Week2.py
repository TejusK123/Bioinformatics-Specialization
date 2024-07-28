import numpy as np 
from collections import deque

'''
Suffix Array Construction Problem: Construct the suffix array of a string.

Input: A string Text.
Output: SuffixArray(Text), as a space-separated collection of integers.
Code Challenge: Solve the Suffix Array Construction Problem.
'''


def suffixarray_construction(text):
	return(' '.join([str(text.index(item)) for item in sorted([text[i:] for i in range(len(text))])]))


'''
Burrows-Wheeler Transform Construction Problem: Construct the Burrows-Wheeler transform of a string.

Input: A string Text.
Output: BWT(Text).
Code Challenge: Solve the Burrows-Wheeler Transform Construction Problem.
'''



def Burrows_Wheeler_Transform(text):
	transform_array = sorted([text[i:] + text[:i]for i in range(len(text))])
	transform = [item[-1] for item in transform_array]
	return(transform)



'''
Inverse Burrows-Wheeler Transform Problem: Reconstruct a string from its Burrows-Wheeler transform.

Input: A string Transform (with a single "
$
$" symbol).
Output: The string Text such that BWT(Text) = Transform.
Code Challenge: Solve the Inverse Burrows-Wheeler Transform Problem.
'''


'''
Helper Function:
Orders instances of nucleotide by position
'''


def order(column1, column2):
	unique, counts = np.unique(column1, return_counts = True)
	hmap = dict(zip(unique, [deque(list(range(item)) + list(range(item))) for item in counts]))
	new_column1 = []
	new_column2 = []
	for i in range(len(column1)):
		new_column1.append(f"{column1[i]}{hmap[column1[i]][0]}")
		hmap[column1[i]].popleft()
	for i in range(len(column2)):
		new_column2.append(f"{column2[i]}{hmap[column2[i]][0]}")
		hmap[column2[i]].popleft()
	return(new_column1, new_column2)


def Inverse_Burrows_Wheelers(text):


	last_array = np.array(list(text))
	first_array = np.array(sorted(text))

	first_column, last_column = order(first_array, last_array)

	testans = ['$0']

	while len(testans) != len(first_column):
		i1 = last_column.index(testans[-1])
		testans.append(first_column[i1])

	ans = [item[0] for item in testans][1:] + ['$']
	return(ans)




'''
Code Challenge: Implement BWMatching.

Input: A string BWT(Text), followed by a space-separated collection of Patterns.
Output: A space-separated list of integers, where the i-th integer corresponds to the number of substring matches of the i-th member of Patterns in Text.
'''


'''
Helper Function
Matches one pattern
'''

def BWMMatch(last_column, Pattern, LastToFirst):
	top = 0
	bottom = len(last_column)
	# print(len(last_column))
	while top <= bottom:
		if len(Pattern) != 0:
			symbol = Pattern[-1]
			Pattern = Pattern[:-1]
			search_range = [item[0] for item in last_column[top:bottom+1]]
			# print(search_range)
			# print(symbol, Pattern)
			if symbol in search_range:
				topIndex = search_range.index(symbol) + top
				bottomIndex = len(search_range) - search_range[::-1].index(symbol) - 1 + top
				# print(bottomIndex, topIndex, bottom, top)
				top = LastToFirst[topIndex]
				bottom = LastToFirst[bottomIndex]
			else:
				return 0
		else:
			return(bottom - top + 1)


#Matches all the patterns
def BWMMatching(Sequence, Patterns):

	last_array = np.array(list(Sequence))
	first_array = np.array(sorted(Sequence))
	first_array, last_array = order(first_array, last_array)
	lasttofirst = [first_array.index(item) for item in last_array]
	ans = []
	for item in Patterns:
		ans.append(str(BWMMatch(last_array, item, lasttofirst)))

	return(ans)


