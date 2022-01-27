"""
problem link: https://rosalind.info/problems/ba7c/
Given: n and a tab-delimited n x n additive matrix.
Return: A weighted adjacency list for the simple tree fitting this matrix.

sample input:
4
0   13  21  22
13  0   12  13
21  12  0   13
22  13  13  0

sample output:
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
"""
def compute_limb_length(node, matrix):
    n = node + 1
    return min(int((matrix[node][i] + matrix[node][j] - matrix[i][j])/2) for i in range(0,n-1) for j in range(0,n-1) if  i != node and j != node and i != j)


class Tree:
    def __init__(self, leave_count):
        self.leave_count = leave_count
        self.tree = {}
        for i in range(0, leave_count):
            self.tree[i] = []

    def add_node(self, father, node, length):
        self.tree[father].append((node, length))
        if node in self.tree:
            self.tree[node].append((father, length))
        else:
            self.tree[node] = [(father, length)]

    def remove_node(self, a, b, length):
        self.tree[a].remove((b, length))
        self.tree[b].remove((a, length))

    def get_path(self, a, b):
        seen = {a:True}
        paths = [[i] for i in self.tree[a]]
        while True:
            path = paths.pop(0)
            last_node_path = path[len(path) - 1][0]
            if last_node_path == b:
                return path
            for edge in self.tree[last_node_path]:
                if edge not in seen:
                    paths.append(path + [edge])
                    seen[edge] = True

    def get_free_node(self):
        return max(self.tree.keys()) + 1

    def print_tree(self):
        for a, tmp in self.tree.items():
            for (b, length) in tmp:
                print(f'{a}->{b}:{length}')


def AdditivePhylogeny(matrix, n, leave_count):
    if n == 2:
        T = Tree(leave_count)
        T.add_node(0, 1, matrix[0][1])
        return T
    limbLength = compute_limb_length(n - 1, matrix)
    for i in range(0, n - 1):
        matrix[i][n - 1] -= limbLength
        matrix[n - 1][i] -= limbLength

    # finding (i,j,n) such that Di,k = Di,n + Dn,k
    i, j = -1, -1
    for k in range(0, n - 1):
        for t in range(0, n - 1):
            if matrix[k][t] == matrix[k][n - 1] + matrix[n - 1][t]:
                i = k
                j = t
                break
        if i != -1:
            break

    x = matrix[i][n - 1]
    T = AdditivePhylogeny(matrix, n - 1, leave_count)

    path = T.get_path(i, j)  # [(node1,len1),(node2,len2),...(j,lenj)]
    length_till_now = 0
    previous_node = i
    free_node = T.get_free_node()
    for (edge, length) in path:
        length_till_now += length
        if length_till_now > x:
            T.remove_node(previous_node, edge, length)
            T.add_node(previous_node, free_node, length-(length_till_now-x))
            T.add_node(free_node, edge, length_till_now - x)
            T.add_node(free_node, n - 1, limbLength)
            break
        if length_till_now == x:
            T.add_node(edge, free_node, x)
            break
        previous_node = edge
    return T


file = open('../Motif Finding/input.txt', 'r')
Lines = file.readlines()
n = int(Lines[0])
matrix = []
for row in Lines[1:]:
    matrix.append(list(map(int, row.strip().split())))
AdditivePhylogeny(matrix, n, n).print_tree()
