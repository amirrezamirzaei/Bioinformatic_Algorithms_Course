"""
problem link: https://rosalind.info/problems/ba5l/
Given: Two long amino acid strings (of length approximately 10,000).
Return: The maximum alignment score of these strings, followed by an alignment achieving this maximum score.
Use the BLOSUM62 scoring matrix and indel penalty σ = 5.

sample input:
PLEASANTLY
MEANLY

sample output:
8
PLEASANTLY
-MEA--N-LY
"""
from utils import blosum62


def align_dynamic(x, y):
    matrix = [[0] * (len(y) + 1) for i in range(len(x) + 1)]
    for i in range(len(y) + 1):
        matrix[0][i] = (i * -5, 'l')
    for i in range(len(x) + 1):
        matrix[i][0] = (i * -5, 'u')
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            m = max(matrix[i][j - 1][0] - 5,
                    matrix[i - 1][j][0] - 5,
                    matrix[i - 1][j - 1][0] + blosum62(y[j - 1], x[i - 1]))
            if m == matrix[i][j - 1][0] - 5:
                matrix[i][j] = (m, 'l')
            elif m == matrix[i - 1][j][0] - 5:
                matrix[i][j] = (m, 'u')
            else:
                matrix[i][j] = (m, 'm')
    aligned_x = ''
    aligned_y = ''
    px = len(x)
    py = len(y)
    while px != 0 or py != 0:
        state = matrix[px][py][1]
        if state == 'm':
            px -= 1
            py -= 1
            aligned_y = y[py] + aligned_y
            aligned_x = x[px] + aligned_x
        elif state == 'l':
            py -= 1
            aligned_y = y[py] + aligned_y
            aligned_x = '-' + aligned_x
        elif state == 'u':
            px -= 1
            aligned_y = '-' + aligned_y
            aligned_x = x[px] + aligned_x
    return aligned_x, aligned_y


def score_last_row_linear_space(x, y):
    previous_row = [i * -5 for i in range(len(y) + 1)]
    for i in range(1, len(x) + 1):
        new_row = [0] * (len(y) + 1)

        for j in range(0, len(y) + 1):

            if j == 0:
                new_row[j] = previous_row[j] - 5
            else:
                m = max(previous_row[j] - 5,
                        new_row[j - 1] - 5,
                        previous_row[j - 1] + blosum62(y[j - 1], x[i - 1]))
                new_row[j] = m
        previous_row = new_row
    return previous_row


def Hirschberg(x, y):  # from http://www.csse.monash.edu.au/~lloyd/tildeAlgDS/Dynamic/Hirsch/
    z = ''
    w = ''
    if len(x) == 0:
        z = z + '-' * len(y)
        w = y
    elif len(y) == 0:
        w = w + '-' * len(x)
        z = x
    elif len(x) == 1 or len(y) == 1:
        z, w = align_dynamic(x, y)
    else:
        xmid = int(len(x) / 2)
        scoreL = score_last_row_linear_space(x[0:xmid], y)
        scoreR = score_last_row_linear_space(x[xmid:][::-1], y[::-1])[::-1]
        t = [sum(x) for x in zip(scoreL, scoreR)]
        ymid = t.index(max(t))
        h1 = Hirschberg(x[0:xmid], y[0:ymid])
        h2 = Hirschberg(x[xmid:], y[ymid:])
        z = h1[0] + h2[0]
        w = h1[1] + h2[1]
    return z, w


with open('input.txt') as f:
    lines = f.readlines()
    v, w = lines[0].strip(), lines[1].strip()

a, b = Hirschberg(v, w)
score = 0
for i, j in zip(a, b):
    if i == '-' or j == '-':
        score -= 5
    else:
        score += blosum62(i, j)
print(score)
print(a)
print(b)
