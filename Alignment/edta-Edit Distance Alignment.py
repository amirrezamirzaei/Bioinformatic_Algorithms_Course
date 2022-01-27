"""
problem link: https://rosalind.info/problems/edta/
Given: Two protein strings s and t in FASTA format (with each string having length at most 1000 aa).
Return: The edit distance dE(s,t) followed by two augmented strings s′ and t′ representing an optimal alignment of s and t.

sample input:
>Rosalind_43
PRETTY
>Rosalind_97
PRTTEIN

sample output:
4
PRETTY--
PR-TTEIN
"""
from utils import read_fasta

s, t = read_fasta('input.txt')
matrix = [[0] * (len(s) + 1) for i in range(len(t) + 1)]

for i in range(len(s) + 1):
    matrix[0][i] = (i, 'l')

for i in range(len(t) + 1):
    matrix[i][0] = (i, 'u')

for i in range(1, len(t) + 1):
    for j in range(1, len(s) + 1):
        m = min(matrix[i - 1][j][0]+1, matrix[i][j - 1][0]+1, matrix[i - 1][j - 1][0] + int(s[j - 1] != t[i - 1]))
        if m == matrix[i - 1][j][0]+1:
            matrix[i][j] = (matrix[i - 1][j][0]+1, 'u')
        elif m == matrix[i - 1][j - 1][0] + int(s[j - 1] != t[i - 1]):
            matrix[i][j] = (matrix[i - 1][j - 1][0] + int(s[j - 1] != t[i - 1]), 'm')
        else:
            matrix[i][j] = (matrix[i][j - 1][0]+1, 'l')


aligned_s = ''
aligned_t = ''
y = len(s)
x = len(t)
while x != 0 or y != 0:
    if matrix[x][y][1] == 'm':
        x -= 1
        y -= 1
        aligned_t = t[x] + aligned_t
        aligned_s = s[y] + aligned_s
    elif matrix[x][y][1] == 'u':
        x -= 1
        aligned_t = t[x] + aligned_t
        aligned_s = '-' + aligned_s
    elif matrix[x][y][1] == 'l':
        y -= 1
        aligned_t = '-' + aligned_t
        aligned_s = s[y] + aligned_s

print(matrix[len(t)][len(s)][0])
print(aligned_s)
print(aligned_t)