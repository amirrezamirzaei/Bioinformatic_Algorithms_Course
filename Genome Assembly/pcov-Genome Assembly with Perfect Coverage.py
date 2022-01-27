"""
problem link: https://rosalind.info/problems/pcov/
Given: A collection of (error-free) DNA k-mers (k≤50) taken from the same strand of a circular chromosome. In this dataset, all k-mers from this strand of the chromosome
are present, and their de Bruijn graph consists of exactly one simple cycle.
Return: A cyclic superstring of minimal length containing the reads (thus corresponding to a candidate cyclic chromosome)

sample input:
ATTAC
TACAG
GATTA
ACAGA
CAGAT
TTACA
AGATT

sample output:
GATTACA
"""
def is_in_circular(s1, s2):
    for i in range(0, len(s2) + 1):
        if s1.startswith(s2[i:]) and s1.endswith(s2[0:i]):
            return True


file = open('../Motif Finding/input.txt', 'r')
reads = list(map(lambda s: s.strip(), file.readlines()))
k = len(reads[0])

graph = {}
k = len(reads[0]) - 1

for read in reads:
    s1 = read[0:k]
    s2 = read[1:k + 1]
    if s1 in graph and s2 not in graph[s1]:
        graph[s1].append(s2)
    elif s1 not in graph:
        graph[s1] = [s2]

assembled_dna = ''
start_node = reads[0][:-1]

while graph:
    assembled_dna += start_node[0]
    v = graph[start_node][0]
    graph[start_node].remove(v)
    if not graph[start_node]:
        del graph[start_node]
    start_node = v
print(assembled_dna)



