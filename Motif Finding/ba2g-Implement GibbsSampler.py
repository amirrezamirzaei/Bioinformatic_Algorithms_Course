"""
problem link: https://rosalind.info/problems/ba2g/
Given: A collection of at most 1000 DNA strings (whose combined length does not exceed 50 kbp).
Return: N50 and N75 for this collection of strings.

sample input:
8 5 100
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
import random


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


def most_probable_kmer_for_motif(dna, profile, pseudo_count=0.1):
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


def gibbs_sampler(Dnas, k, t, N=1000, random_start=2000):
    best_motifs = None
    best_score = 0
    for iter in range(0, random_start):
        motifs = [s[random.randrange(len(s) - k + 1):][0:k] for s in Dnas]
        best_iter_motif = motifs
        best_iter_score = motif_score(motifs)
        for j in range(0, N):
            i = random.randrange(t)
            profile = make_profile(motifs[0:i] + motifs[i + 1:])
            motifs[i] = most_probable_kmer_for_motif(Dnas[i], profile, pseudo_count=0.5)
            score = motif_score(motifs)
            if score > best_iter_score:
                best_iter_motif = motifs
                best_iter_score = score
        if best_iter_score > best_score:
            best_score = best_iter_score
            best_motifs = best_iter_motif
    return best_motifs


file = open('input.txt', 'r')
Lines = file.readlines()
k, t, N = list(map(int, Lines[0].split()))
Dna = list(map(lambda s: s.strip(), Lines[1:]))
for motif in gibbs_sampler(Dna, k, t, N=N, random_start=20):
    print(motif)
