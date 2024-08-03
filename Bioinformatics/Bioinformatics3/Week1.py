import numpy as np 
import networkx as nx


'''
Code Challenge: Solve the Change Problem. The DPChange pseudocode is reproduced below for your convenience.

Input: An integer money and an array Coins = (coin1, ..., coind).
Output: The minimum number of coins with denominations Coins that changes money.
'''

def DPChange(money, Coins):
	Coins = sorted(Coins, reverse = True)
	MinNumCoins = [0] * (money + 1)  
	MinNumCoins[0] = 0
	for m in range(1,money + 1):
		MinNumCoins[m] = float('inf')
		for i in range(len(Coins)):
			if m >= Coins[i]:
				if MinNumCoins[m - Coins[i]] + 1 < MinNumCoins[m]:
					MinNumCoins[m] = MinNumCoins[m - Coins[i]] + 1
	return(MinNumCoins[-1])




'''
Code Challenge: Find the length of a longest path in the Manhattan Tourist Problem.

Input: Integers n and m, followed by an n × (m + 1) matrix Down and an (n + 1) × m matrix Right. The two matrices are separated by the "-" symbol.
Output: The length of a longest path from source (0, 0) to sink (n, m) in the rectangular grid whose edges are defined by the matrices Down and Right.
'''

def ManhattenTourist(n, m, Down, Right, Diag = []):
	s = np.zeros((Right.shape[0], Down.shape[1]))
	s[0,0] = 0
	for i in range(1, n + 1):
		s[i,0] = s[i-1,0] + Down[i-1,0]
	for j in range(1, m + 1):
		s[0,j] = s[0, j-1] + Right[0, j-1]
	
	if len(Diag) == 0:
		for i in range(1,n+1):
			for j in range(1, m+1):
				s[i,j] = max(s[i-1,j] + Down[i-1,j], s[i,j-1] + Right[i,j-1])
	else:
		for i in range(1,n+1):
			for j in range(1, m+1):
				s[i,j] = max(max(s[i-1,j] + Down[i-1,j], s[i,j-1] + Right[i,j-1]), s[i-1,j-1] + Diag[i-1,j-1])


	return(s)



#Usage

Down = np.array([[1,0,2,4,3],[4,6,5,2,1],[4,4,5,2,1],[5,6,8,5,3]])
Right = np.array([[3,2,4,0],[3,2,4,2],[0,7,3,3],[3,3,0,2],[1,3,2,2]])

print(ManhattenTourist(4,4,Down,Right)[-1,-1])


'''
Code Challenge: Use OutputLCS (reproduced below) to solve the Longest Common Subsequence Problem.

Input: Two strings s and t.
Output: A longest common subsequence of s and t. (Note: more than one solution may exist, in which case you may output any one.)
OutputLCS(backtrack, v, i, j)
    if i = 0 or j = 0
        return ""
    if backtracki, j = "↓"
        return OutputLCS(backtrack, v, i - 1, j)
    else if backtracki, j = "→"
        return OutputLCS(backtrack, v, i, j - 1)
    else
        return OutputLCS(backtrack, v, i - 1, j - 1) + vi
'''

'''
Helper Functions:
LCSBackTrack:
returns the backtrack matrix so that you can calculate the Longest Common Subsequence
'''
def LCSBackTrack(v, w):
	n = len(v)
	m = len(w)
	s = np.zeros((n+1,m+1))
	backtrack = np.zeros((n+1,m+1))
	for i in range(1, n+1):
		for j in range(1, m+1):
			match = 0 
			if v[i-1] == w[j-1]:
				match = 1 
			s[i,j] = max(s[i-1,j], s[i,j-1], s[i-1,j-1] + match)
			if s[i, j] == s[i-1, j]:
				backtrack[i, j] = 1 #down
			elif s[i, j] == s[i, j-1]:
				backtrack[i, j] = 2 #right
			elif s[i, j] == s[i-1, j-1] + 1 and v[i-1] == w[j-1]:
				backtrack[i, j] = 3 #Diag
	return backtrack

def OutputLCS(v, w):
	backtrack = LCSBackTrack(v,w)
	i = len(v)
	j = len(w)
	lcs = ''
	while i > 0 and j > 0:
		if 1 == backtrack[i, j]:
			i -= 1
		elif 2 == backtrack[i, j]:
			j -= 1
		else:
			i -= 1
			j -= 1
			lcs += v[i]

	return(lcs[::-1])



print(OutputLCS('AACCTTGG', 'ACACTGTGA'))



'''

Code Challenge: Solve the Longest Path in a DAG Problem.

Input: An integer representing the starting node to consider in a graph, followed by an integer representing the ending node to consider, followed by a list of edges in the graph. The edge notation "0 1 7" indicates that an edge connects node 0 to node 1 with weight 7.  You may assume a given topological order corresponding to nodes in increasing order.
Output: The length of a longest path in the graph, followed by a longest path as a sequence of space-separated node labels. (If multiple longest paths exist, you may return any one.)

'''

def longest_path(G, s, sink):
    
    start_node = s
    topoorder = nx.topological_sort(G)
    dist = dict.fromkeys(G.nodes, -float('inf'))
    dist[s] = 0
    for n in topoorder:
        for s in G.successors(n):
            
            if dist[s] < dist[n] + G.edges[n,s]['weight']:

                dist[s] = dist[n] + G.edges[n,s]['weight']

    
    path = [sink]
    ssink = sink 
    while start_node not in (G.predecessors(ssink)):
        cur_max = 0
        
        for item in G.predecessors(ssink):
            if dist[item] > cur_max:
                cur_max = dist[item]
                new = item
        
        path.append(new)
        ssink = new

    path.append(start_node)

    #My manual implementation works like 95% of the time I don't know what the edge case is
    
    return(dist[sink], nx.dag_longest_path(G))




'''
Usage
Sample Input:
0 4
0 1 7
0 2 4
2 3 2
1 4 1
3 4 3
'''


#put your own testset
with open(r"LongestPathtest.txt") as f:
    lines = f.readlines()

n, m = int(lines[0].split(' ')[0]), int(lines[0].split(' ')[1][:-1])
print(m)
G = nx.DiGraph()
for line in lines[1:]:
    items = line.split(' ')
    items[-1] = items[-1][:-1] if items[-1][-1] == '\n' else items[-1]
    items = [int(item) for item in items]
    # print(items)


    if [items[0],items[1]] in G.edges:
        if items[2] > G.edges[items[0],items[1]]['weight']: # in the case of duplicate edges, pick the edge with the largest weight to add to the graph
            G.remove_edge(items[0],items[1])
            G.add_edge(items[0],items[1],weight = items[2])
    else:
        G.add_edge(items[0], items[1], weight = items[2])


print(longest_path(G,n,m))
