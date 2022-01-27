"""
problem link: https://rosalind.info/problems/sims/
Given: Two DNA strings s and t, where s has length at most 10
kbp and t represents a motif of length at most 1 kbp.
Return: An optimal fitting alignment score with respect to the
mismatch score defined above, followed by an optimal fitting
alignment of a substring of s against t. If multiple
such alignments exist, then you may output any one.

sample input:
>Rosalind_54
GCAAACCATAAGCCCTACGTGCCGCCTGTTTAAACTCGCGAACTGAATCTTCTGCTTCACGGTGAAAGTACCACAATGGTATCACACCCCAAGGAAAC
>Rosalind_46
GCCGTCAGGCTGGTGTCCG

sample output:
5
ACCATAAGCCCTACGTG-CCG
GCCGTCAGGC-TG-GTGTCCG
"""
from utils import read_fasta

s, t = read_fasta('input.txt')

matrix = [[0] * (len(s) + 1) for i in range(len(t) + 1)]

for i in range(len(s) + 1):
    matrix[0][i] = (-i, 'u')

for i in range(len(t) + 1):
    matrix[i][0] = (-i, 'l')

for i in range(1, len(t) + 1):
    for j in range(1, len(s) + 1):
        match_score = matrix[i - 1][j - 1][0] + (1 if int(s[j - 1] == t[i - 1]) else -1)
        m = max(matrix[i - 1][j][0] - 1, matrix[i][j - 1][0] - 1,
                match_score, (1 if int(s[j - 1] == t[i - 1]) else -1) - (i - 1))
        if m == match_score:
            matrix[i][j] = (m, 'm')
        elif m == matrix[i - 1][j][0] - 1:
            matrix[i][j] = (m, 'u')
        elif m == matrix[i][j - 1][0] - 1:
            matrix[i][j] = (m, 'l')
        elif m == (1 if int(s[j - 1] == t[i - 1]) else -1) - (i - 1):
            matrix[i][j] = (m, 'b')

alignment_max, x_max, y_max = 0, len(t), 0
for i in range(len(s)):
    if matrix[x_max][i][0] > alignment_max:
        alignment_max = matrix[x_max][i][0]
        y_max = i

aligned_s = ''
aligned_t = ''
flag = True
while (x_max != 0 or y_max != 0) and flag:
    flag = (matrix[x_max][y_max][1] != 'b')
    if matrix[x_max][y_max][1] == 'm' or matrix[x_max][y_max][1] == 'b':
        x_max -= 1
        y_max -= 1
        aligned_t = t[x_max] + aligned_t
        aligned_s = s[y_max] + aligned_s
    elif matrix[x_max][y_max][1] == 'u':
        x_max -= 1
        aligned_t = t[x_max] + aligned_t
        aligned_s = '-' + aligned_s
    elif matrix[x_max][y_max][1] == 'l':
        y_max -= 1
        aligned_t = '-' + aligned_t
        aligned_s = s[y_max] + aligned_s
aligned_t = t[0:x_max] + aligned_t

print(alignment_max)
print(aligned_s)
print(aligned_t)
