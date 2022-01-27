"""
problem link: https://rosalind.info/problems/laff/
Given: Two protein strings s and t in FASTA format (each having length at most 10,000 aa).
Return: The maximum local alignment score of s and t, followed by substrings r and u of s and t, respectively, that correspond to the optimal local alignment of s and t. Use:
    The BLOSUM62 scoring matrix.
    Gap opening penalty equal to 11.
    Gap extension penalty equal to 1.

sample input:
>Rosalind_8
PLEASANTLY
>Rosalind_18
MEANLY

sample output:
12
LEAS
MEAN
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
best_score, best_x, best_y = 0, 0, 0

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
                matrix_down[i][j][0], blosum62(s[j - 1], t[i - 1]))

        if m == matrix[i - 1][j - 1][0] + blosum62(s[j - 1], t[i - 1]):
            matrix[i][j] = (m, 'm')
        elif m == matrix_right[i][j][0]:
            matrix[i][j] = (m, 'l')
        elif m == matrix_down[i][j][0]:
            matrix[i][j] = (m, 'u')
        elif m == blosum62(s[j - 1], t[i - 1]):
            matrix[i][j] = (m, 'b')

        if m > best_score:
            best_score = m
            best_x = i
            best_y = j

aligned_s = ''
aligned_t = ''

state = 'm'
flag = True
while (best_x != 0 or best_y != 0) and flag:
    flag = (state != 'b')
    if state == 'm' or state == 'b':

        best_x -= 1
        best_y -= 1
        aligned_t = t[best_x] + aligned_t
        aligned_s = s[best_y] + aligned_s
        state = matrix[best_x][best_y][1]

    elif state == 'l':
        state = matrix_right[best_x][best_y][1]
        best_y -= 1
        aligned_s = s[best_y] + aligned_s

    elif state == 'u':
        state = matrix_down[best_x][best_y][1]
        best_x -= 1
        aligned_t = t[best_x] + aligned_t

print(best_score)
print(aligned_s)
print(aligned_t)
