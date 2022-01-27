"""
problem link: https://rosalind.info/problems/ba7a/
Given: An integer n followed by the adjacency list of a weighted tree with n leaves.
Return: A space-separated n x n (di, j), where di, j is the length of the path between leaves i and j.

sample input:
4
0->4:11
1->4:2
2->5:6
3->5:7
4->0:11
4->1:2
4->5:4
5->4:4
5->3:7
5->2:6

sample output:
0   13  21  22
13  0   12  13
21  12  0   13
22  13  13  0
"""
file1 = open('../Motif Finding/input.txt', 'r')
Lines = file1.readlines()
n = int(Lines[0])
tree = {}
for edge in Lines[1:]:
    edge = edge.strip().split(':')
    nodes = edge[0].split('->')
    a = int(nodes[0])
    b = int(nodes[1])
    if b < a:
        continue
    length = int(edge[1])
    if b in tree:
        tree[b].append((a, length))
    else:
        tree[b] = [(a, length)]
    if a in tree:
        tree[a].append((b, length))
    else:
        tree[a] = [(b, length)]

output_matrix = [[0 for j in range(0, n)] for i in range(0, n)]
for i in range(0, n):
    queue = [(i, 0)]
    seen = {i: True}
    while queue:
        (edge, length) = queue.pop()
        if edge < n and edge != i:
            output_matrix[i][edge] = length
        else:
            for edge_tmp, length_tmp in tree[edge]:
                if edge_tmp not in seen:
                    seen[edge_tmp] = True
                    queue.append((edge_tmp, length_tmp + length))

for l in output_matrix:
    print(str(l).replace('[', '').replace(']', '').replace(',', ''))
