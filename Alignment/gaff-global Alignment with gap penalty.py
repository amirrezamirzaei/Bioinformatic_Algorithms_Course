"""
problem link: https://rosalind.info/problems/gaff/
Given: Two protein strings s and t in FASTA format (each of length at most 100 aa).
Return: The maximum alignment score between s and t, followed by two augmented strings s′ and t′ representing an optimal alignment of s and t. Use:
    The BLOSUM62 scoring matrix.
    Gap opening penalty equal to 11.
    Gap extension penalty equal to 1.

sample input:
>Rosalind_49
PRTEINS
>Rosalind_47
PRTWPSEIN

sample output:
8
PRT---EINS
PRTWPSEIN-
"""
from utils import read_fasta, blosum62

s, t = read_fasta('input.txt')

gap_opening = -11
gap_extension = -1
matrix = [[0] * (len(s) + 1) for i in range(len(t) + 1)]
matrix_right = [[0] * (len(s) + 1) for i in range(len(t) + 1)]
matrix_down = [[0] * (len(s) + 1) for i in range(len(t) + 1)]

matrix[0][0] = (0, 'u')
matrix_down[0][0] = (0, 'u')
matrix_right[0][0] = (0, 'l')

for i in range(1, len(t) + 1):
    matrix[i][0] = (-11 - (i - 1), 'u')
    matrix_right[i][0] = (-9999, 'l')
    matrix_down[i][0] = (-9999, 'u')
for j in range(1, len(s) + 1):
    matrix[0][j] = (-11 - (j - 1), 'l')
    matrix_right[0][j] = (-9999, 'l')
    matrix_down[0][j] = (-9999, 'u')

for i in range(1, len(t) + 1):
    for j in range(1, len(s) + 1):
        max_right = max(gap_extension + matrix_right[i][j - 1][0], gap_opening + matrix[i][j - 1][0])
        if max_right == gap_opening + matrix[i][j - 1][0]:
            matrix_right[i][j] = (gap_opening + matrix[i][j - 1][0], 'm')
        else:
            matrix_right[i][j] = (gap_extension + matrix_right[i][j - 1][0], 'l')

        max_left = max(gap_extension + matrix_down[i - 1][j][0], gap_opening + matrix[i - 1][j][0])
        if max_left == gap_opening + matrix[i - 1][j][0]:
            matrix_down[i][j] = (gap_opening + matrix[i - 1][j][0], 'm')
        else:
            matrix_down[i][j] = (gap_extension + matrix_down[i - 1][j][0], 'u')

        m = max(matrix[i - 1][j - 1][0] + blosum62(s[j - 1], t[i - 1]), matrix_right[i][j][0],
                matrix_down[i][j][0])

        if m == matrix[i - 1][j - 1][0] + blosum62(s[j - 1], t[i - 1]):
            matrix[i][j] = (m, 'm')
        elif m == matrix_right[i][j][0]:
            matrix[i][j] = (m, 'l')
        elif m == matrix_down[i][j][0]:
            matrix[i][j] = (m, 'u')

aligned_s = ''
aligned_t = ''
y = len(s)
x = len(t)
alignment_score = max(matrix[x][y][0], matrix_right[x][y][0], matrix_down[x][y][0])
state = 'm' if alignment_score == matrix[x][y][0] else ('l' if alignment_score == matrix_right[x][y][0] else 'u')
while x != 0 or y != 0:
    if state == 'm':

        x -= 1
        y -= 1
        aligned_t = t[x] + aligned_t
        aligned_s = s[y] + aligned_s
        state = matrix[x][y][1]

    elif state == 'l':
        state = matrix_right[x][y][1]
        y -= 1
        aligned_t = '-' + aligned_t
        aligned_s = s[y] + aligned_s

    elif state == 'u':
        state = matrix_down[x][y][1]
        x -= 1
        aligned_t = t[x] + aligned_t
        aligned_s = '-' + aligned_s

print(matrix[len(t)][len(s)][0])
print(aligned_s)
print(aligned_t)
