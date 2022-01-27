"""
problem link: https://rosalind.info/problems/ba2e/
Given: Integers k and t, followed by a collection of strings Dna.
Return: A collection of strings BestMotifs resulting from running GreedyMotifSearch(Dna, k, t) with pseudocounts.
If at any step you find more than one Profile-most probable k-mer in a given string, use the one occurring first.

sample input:
3 5
GGCGTTCAGGCA
AAGAATCAGTCA
CAAGGAGTTCGC
CACGTCAATCAC
CAATAATATTCG

sample output:
TTC
ATC
TTC
ATC
TTC
"""


def motif_score(motifs):
    score = 0
    for i in range(0, len(motifs[0])):
        count = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        for motif in motifs:
            count[motif[i]] += 1
        score += count[max(count, key=count.get)]
    return score


def make_profile(motifs):
    profile = []
    for i in range(0, len(motifs[0])):
        count = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        for motif in motifs:
            count[motif[i]] += 1
        profile.append(count)
    return profile


def most_probable_kmer_for_motif(dna, profile, pseudo_count):
    k = len(profile)
    best_kmer = dna[0:k]
    max_prob = 0
    for i in range(0, len(dna) - k + 1):
        kmer = dna[i:i + k]
        prob = 1
        for i in range(0, k):
            prob *= (profile[i][kmer[i]] + pseudo_count)
        if prob > max_prob:
            max_prob = prob
            best_kmer = kmer

    return best_kmer


def greedy_motif_search(Dnas, k, t, pseudo_count=1):
    best_score = 0
    best_motifs = [s[0:k] for s in Dnas]

    for j in range(0, len(Dnas[0]) - k + 1):
        motifs = [Dnas[0][j:j + k]]
        for i in range(1, t):
            profile = make_profile(motifs)
            motif_i = most_probable_kmer_for_motif(Dnas[i], profile, pseudo_count)
            motifs.append(motif_i)
        score = motif_score(motifs)
        if score > best_score:
            best_score = score
            best_motifs = motifs
    return best_motifs


file = open('input.txt', 'r')
Lines = file.readlines()
k, t = list(map(int, Lines[0].split()))
Dna = list(map(lambda s: s.strip(), Lines[1:]))
for motif in greedy_motif_search(Dna, k, t):
    print(motif)
