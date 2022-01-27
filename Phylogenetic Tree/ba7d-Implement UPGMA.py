"""
problem link: https://rosalind.info/problems/ba7d/
Given: An integer n followed by a space-delimited n x n distance matrix.
Return: An adjacency list for the ultrametric tree output by UPGMA. Weights should be accurate to three decimal places.

sample input:
4
0   20  17  11
20  0   20  13
17  20  0   10
11  13  10  0

sample output:
0->5:7.000
1->6:8.833
2->4:5.000
3->4:5.000
4->2:5.000
4->3:5.000
4->5:2.000
5->0:7.000
5->4:2.000
5->6:1.833
6->5:1.833
6->1:8.833
"""
import numpy as np


class Tree:
    def __init__(self, leave_count):
        self.leave_count = leave_count
        self.tree = {}
        self.node_age = {}
        for i in range(0, leave_count):
            self.tree[i] = []
            self.node_age[i] = 0
        self.last_node_num = leave_count - 1

    def add_node(self, father, child, age):
        self.node_age[father] = age
        if father in self.tree:
            self.tree[father].append(child)
        else:
            self.tree[father] = [child]
        self.tree[child].append(father)

    def get_free_node(self):
        self.last_node_num += 1
        return self.last_node_num

    def print_tree(self):
        for a, tmp in self.tree.items():
            for b in tmp:
                print(f'{a}->{b}:{np.abs((self.node_age[b] - self.node_age[a])):.3f}')


def UPGMA(matrix, n):
    clusters = [[i] for i in range(0, n)]
    matrix_rows = [i for i in range(0, n)]
    T = Tree(n)
    original_mat = matrix.copy()
    while len(matrix_rows) != 1:
        # finding the closest clusters
        minimum, i, j = None, None, None
        for x in range(0, len(matrix_rows)):
            for y in range(0, len(matrix_rows)):
                if (x != y) and (not minimum or matrix[x][y] < minimum):
                    minimum = matrix[x][y]
                    i = x
                    j = y

        m = T.get_free_node()
        age_m = matrix[i][j] / 2.0
        new_cluster = clusters[i] + clusters[j]
        del clusters[i]
        del clusters[j - 1 if j > i else j]
        clusters.append(new_cluster)
        T.add_node(m, matrix_rows[i], age_m)
        T.add_node(m, matrix_rows[j], age_m)
        del matrix_rows[i]
        del matrix_rows[j - 1 if j > i else j]
        matrix_rows.append(m)
        new_row = np.zeros(len(matrix_rows))
        index = 0
        for t in range(0, len(matrix_rows)+1):
            if t != i and t != j:
                for a in clusters[index]:
                    for b in clusters[-1]:
                        new_row[index] += original_mat[a][b]
                new_row[index] /= (len(clusters[index]) * len(clusters[-1]))
                index += 1
        matrix = np.delete(matrix, (i, j), 0)
        matrix = np.delete(matrix, (i, j), 1)
        matrix = np.vstack((matrix, new_row[0:-1]))
        matrix = np.concatenate((matrix, new_row.reshape(-1, 1)), axis=1)
    return T


file = open('../Motif Finding/input.txt', 'r')
Lines = file.readlines()
n = int(Lines[0])
matrix = np.zeros((n, n))
for index, row in enumerate(Lines[1:]):
    matrix[index] = list(map(int, row.strip().split()))

UPGMA(matrix, n).print_tree()
