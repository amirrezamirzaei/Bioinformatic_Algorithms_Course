"""
problem link: https://rosalind.info/problems/ba2f/
Given: Positive integers k and t, followed by a collection of strings Dna.
Return: A collection BestMotifs resulting from running RandomizedMotifSearch(Dna, k, t) 1000 times.
Remember to use pseudocounts!

sample input:
8 5
CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA
GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
TAGTACCGAGACCGAAAGAAGTATACAGGCGT
TAGATCAAGTTTCAGGTGCACGTCGGTGAACC
AATCCACCAGCTCCACGTGCAATGTTGGCCTA

sample output:
TCTCGGGG
CCAAGGTG
TACAGGCG
TTCAGGTG
TCCACGTG
"""
from random import randrange


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


def random_motif_search(Dnas, k, t, pseudo_count=1, iteration=1):
    best_motif = None
    best_score = 0
    dna_length = len(Dnas[0])
    for iter in range(0, iteration):
        motifs = [s[randrange(dna_length - k + 1):][0:k] for s in Dnas]
        best_iter_motif = motifs
        best_iter_score = 0
        while True:
            profile = make_profile(motifs)
            motifs = [most_probable_kmer_for_motif(dna, profile, pseudo_count) for dna in Dnas]
            score = motif_score(motifs)
            if score > best_iter_score:
                best_iter_score = score
                best_iter_motif = motifs
            else:
                if best_iter_score > best_score:
                    best_score = best_iter_score
                    best_motif = best_iter_motif
                break
    return best_motif


file = open('input.txt', 'r')
Lines = file.readlines()
k, t = list(map(int, Lines[0].split()))
Dna = list(map(lambda s: s.strip(), Lines[1:]))
for motif in random_motif_search(Dna, k, t, pseudo_count=1, iteration=1000):
    print(motif)
