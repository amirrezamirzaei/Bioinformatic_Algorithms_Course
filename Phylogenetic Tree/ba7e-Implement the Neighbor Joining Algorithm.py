"""
problem link: https://rosalind.info/problems/ba7e/
Given: An integer n, followed by a space-separated n x n distance matrix.
Return: An adjacency list for the tree resulting from applying the neighbor-joining algorithm.
Edge-weights should be accurate to two decimal places (they are provided to three decimal places in the sample output below).

sample input:
4
0   23  27  20
23  0   30  28
27  30  0   30
20  28  30  0

sample output:
0->4:8.000
1->5:13.500
2->5:16.500
3->4:12.000
4->5:2.000
4->0:8.000
4->3:12.000
5->1:13.500
5->2:16.500
5->4:2.000
"""
import numpy as np


class Tree:
    def __init__(self, leave_count):
        self.leave_count = leave_count
        self.tree = {}
        for i in range(0, leave_count):
            self.tree[i] = []
        self.last_node_num = leave_count - 1

    def add_node(self, father, node, length):
        if father in self.tree:
            self.tree[father].append((node, length))
        else:
            self.tree[father] = [(node, length)]
        if node in self.tree:
            self.tree[node].append((father, length))
        else:
            self.tree[node] = [(father, length)]

    def get_free_node(self):
        self.last_node_num += 1
        return self.last_node_num

    def print_tree(self):
        for a, tmp in self.tree.items():
            for (b, length) in tmp:
                print(f'{a}->{b}:{length:.3f}')


def get_neighbors(matrix):
    n = matrix.shape[0]
    R = np.array([np.sum(matrix[i]) for i in range(0, n)])
    neighbor_joining_matrix = (((n - 2) * matrix).T - R).T - R

    i, j, minimum = None, None, None
    for x in range(0, n):
        for y in range(0, n):
            if (x != y) and ((not minimum) or (neighbor_joining_matrix[x][y] < minimum)):
                minimum = neighbor_joining_matrix[x][y]
                i = x
                j = y
    return i, j, R[i], R[j]


def NeighborJoining(matrix, n, tree, matrix_row_nodes):
    if n == 2:
        tree.add_node(matrix_row_nodes[0], matrix_row_nodes[1], matrix[0][1])
        return tree
    i, j, dis_i, dis_j = get_neighbors(matrix)
    sigma = (dis_i - dis_j) / (n - 2)
    limb_length_i = 0.5 * (matrix[i][j] + sigma)
    limb_length_j = 0.5 * (matrix[i][j] - sigma)
    distance_m = np.zeros((n - 1))
    index = 0
    for k in range(0, n):
        if k != i and k != j:
            distance_m[index] = (matrix[k][i] + matrix[k][j] - matrix[i][j]) / 2
            index += 1
    matrix = np.delete(matrix, (i, j), 0)
    matrix = np.delete(matrix, (i, j), 1)
    matrix = np.vstack((matrix, distance_m[0:-1]))
    matrix = np.concatenate((matrix, distance_m.reshape(-1, 1)), axis=1)

    new_rows = matrix_row_nodes.copy()
    del new_rows[i]
    del new_rows[j - 1 if j > i else j]
    m = tree.get_free_node()
    new_rows.append(m)
    tree = NeighborJoining(matrix, n - 1, tree, new_rows)
    tree.add_node(m, matrix_row_nodes[i], limb_length_i)
    tree.add_node(m, matrix_row_nodes[j], limb_length_j)
    return tree


file = open('../Motif Finding/input.txt', 'r')
Lines = file.readlines()
n = int(Lines[0])
matrix = np.zeros((n, n))
for index, row in enumerate(Lines[1:]):
    matrix[index] = list(map(int, row.strip().split()))

NeighborJoining(matrix, n, Tree(n), [i for i in range(0, n)]).print_tree()
