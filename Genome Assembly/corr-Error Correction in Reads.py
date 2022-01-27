"""
problem link: https://rosalind.info/problems/corr/
Given: A collection of up to 1000 reads of equal length (at most 50 bp) in FASTA format.
Some of these reads were generated with a single-nucleotide error. For each read s in the dataset, one of the following applies:
    s was correctly sequenced and appears in the dataset at least twice (possibly as a reverse complement);
    s is incorrect, it appears in the dataset exactly once, and its Hamming distance is 1 with respect to exactly one correct read in the dataset (or its reverse complement).
Return: A list of all corrections in the form "[old read]->[new read]".
(Each correction must be a single symbol substitution, and you may return the corrections in any order.)

sample input:
>Rosalind_52
TCATC
>Rosalind_44
TTCAT
>Rosalind_68
TCATC
>Rosalind_28
TGAAA
>Rosalind_95
GAGGA
>Rosalind_66
TTTCA
>Rosalind_33
ATCAA
>Rosalind_21
TTGAT
>Rosalind_18
TTTCC

sample output:
TTCAT->TTGAT
GAGGA->GATGA
TTTCC->TTTCA
"""


def get_reverse_complement(s):
    rs = ''
    for i in s:
        if i == 'A':
            rs = 'T' + rs
        elif i == 'T':
            rs = 'A' + rs
        elif i == 'C':
            rs = 'G' + rs
        elif i == 'G':
            rs = 'C' + rs
    return rs


def hamming_distance(s, t):
    distance = 0
    for c1, c2 in zip(s, t):
        if c1 != c2:
            distance += 1
    return distance


dataset = {}

file1 = open('../Motif Finding/input.txt', 'r')
Lines = file1.readlines()
for dna in Lines:
    dna = dna.strip()
    if dna.startswith('>'):
        continue
    if dna in dataset:
        dataset[dna] += 1
    else:
        dataset[dna] = 1

errors = []
for dna, repeat in dataset.items():
    rdna = get_reverse_complement(dna)
    if repeat > 1 or rdna in dataset:
        continue

    for fix, repeat2 in dataset.items():
        if repeat2 > 1 or get_reverse_complement(fix) in dataset:
            if hamming_distance(dna, fix) == 1:
                errors.append((dna, fix))
                break
            elif hamming_distance(rdna, fix) == 1:
                errors.append((dna, get_reverse_complement(fix)))
                break

for (dna, correct) in errors:
    print(f'{dna}->{correct}')
