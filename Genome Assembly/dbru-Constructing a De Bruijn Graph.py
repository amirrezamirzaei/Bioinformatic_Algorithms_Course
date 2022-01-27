"""
problem link: https://rosalind.info/problems/dbru/
Given: A collection of up to 1000 (possibly repeating) DNA strings of equal length (not exceeding 50 bp)
corresponding to a set S of (k+1)-mers.
Return: The adjacency list corresponding to the de Bruijn graph corresponding to S∪Src.

sample input:
TGAT
CATG
TCAT
ATGC
CATC
CATC

sample output:
(ATC, TCA)
(ATG, TGA)
(ATG, TGC)
(CAT, ATC)
(CAT, ATG)
(GAT, ATG)
(GCA, CAT)
(TCA, CAT)
(TGA, GAT)
"""
file = open('../Motif Finding/input.txt', 'r')
reads = list(map(lambda s: s.strip(), file.readlines()))
file.close()
reverse_reads = []
t = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
for dna in reads:
    reverse = ''
    for i in dna:
        reverse = t[i] + reverse
    reverse_reads.append(reverse)
reads = reads + reverse_reads
graph = {}
k = len(reads[0]) - 1
while reads:
    dna = reads.pop()
    s1 = dna[0:k]
    s2 = dna[1:k + 1]
    if s1 in graph and s2 not in graph[s1]:
        graph[s1].append(s2)
    elif s1 not in graph:
        graph[s1] = [s2]

for s1, adjacency_list in sorted(graph.items()):
        for s2 in adjacency_list:
            print(f'({s1}, {s2})')