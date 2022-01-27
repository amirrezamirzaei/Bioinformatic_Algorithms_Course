"""
problem link: https://rosalind.info/problems/ba2b/
Given: An integer k and a collection of strings Dna.
Return: A k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern. (If multiple answers exist, you may return any one.)

sample input:
3
AAATTGACGCAT
GACGACCACGTT
CGTCAGCGCCTG
GCTGAGCACCGG
AGTACGGGACAG

sample output:
GAC
"""

def get_kmer_from_number(n, k):
    binary_num = bin(n)[2:].zfill(2 * k)
    kmer = ''
    for i in range(0, k * 2, 2):
        if binary_num[i:i + 2] == '00':
            kmer += 'A'
        elif binary_num[i:i + 2] == '01':
            kmer += 'T'
        elif binary_num[i:i + 2] == '10':
            kmer += 'C'
        elif binary_num[i:i + 2] == '11':
            kmer += 'G'
    return kmer


def compute_max_hamming_dist(kmer, dna):
    max_hamming = 0
    for i in range(0, len(dna) - len(kmer) + 1):
        tmp = dna[i:i + len(kmer)]
        hamming = 0
        for i in range(0, len(kmer)):
            if kmer[i] == tmp[i]:
                hamming += 1
        if hamming > max_hamming:
            max_hamming = hamming
            if max_hamming == len(kmer):
                break
    return max_hamming


file = open('input.txt', 'r')
Lines = file.readlines()
k = int(Lines[0])
dnas = list(map(lambda s: s.strip(), Lines[1:]))

max_score = 0
max_kmer = ''
for i in range(0, 4 ** k):
    kmer = get_kmer_from_number(i, k)
    score = 0
    for dna in dnas:
        score += compute_max_hamming_dist(kmer, dna)
    if score > max_score:
        max_score = score
        max_kmer = kmer
print(max_kmer)
