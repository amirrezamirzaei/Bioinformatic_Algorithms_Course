"""
problem link: https://rosalind.info/problems/asmq/
Given: A collection of at most 1000 DNA strings (whose combined length does not exceed 50 kbp).
Return: N50 and N75 for this collection of strings.

sample input:
GATTACA
TACTACTAC
ATTGAT
GAAGA

sample output:
7 6
"""


def compute_Nxx_statistic(reads, xx):
    reads_length = list(map(len, reads))
    reads_length.sort(reverse=True)
    all_bp = sum(reads_length)
    for L in reads_length:
        bps = 0
        for read in reads:
            if len(read) >= L:
                bps += len(read)
        if (bps / all_bp) * 100 >= xx:
            return L


file = open('../Motif Finding/input.txt', 'r')
dna = list(map(lambda s: s.strip(), file.readlines()))
print(f'{compute_Nxx_statistic(dna, 50)} {compute_Nxx_statistic(dna, 75)}')
