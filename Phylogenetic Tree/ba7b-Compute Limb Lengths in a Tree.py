"""
problem link: https://rosalind.info/problems/ba7b/
Given: An integer n, followed by an integer j between 0 and n - 1, followed by a space-separated additive distance matrix D (whose elements are integers).
Return: The limb length of the leaf in Tree(D) corresponding to row j of this distance matrix (use 0-based indexing).

sample input:
4
1
0   13  21  22
13  0   12  13
21  12  0   13
22  13  13  0

sample output:
2
"""
file = open('../Motif Finding/input.txt', 'r')
Lines = file.readlines()
n = int(Lines[0])
node = int(Lines[1])
matrix = []
for row in Lines[2:]:
    matrix.append(list(map(int, row.strip().split())))

limb_len = min(int((matrix[node][i] + matrix[node][j] - matrix[i][j])/2) for i in range(0,n-1) for j in range(0,n-1) if  i != node and j != node and i != j)
print(limb_len)
