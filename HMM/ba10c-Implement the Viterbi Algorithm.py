"""
problem link: https://rosalind.info/problems/ba10c/
Given: A string x, followed by the alphabet Σ from which x was constructed, followed by the states States, transition matrix Transition, and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
Return: A path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π.

sample input:
xyxzzxyxyy
--------
x   y   z
--------
A   B
--------
    A   B
A   0.641   0.359
B   0.729   0.271
--------
    x   y   z
A   0.117   0.691   0.192
B   0.097   0.42    0.483

sample output:
AAABBAAAAA
"""
from io import StringIO
import numpy as np
import pandas as pd

class HMM:
    def __init__(self, emission, transition, alphabet, states):
        self.emission = emission
        self.transition = transition
        self.alphabet = alphabet
        self.states = states

    def viterbi(self, evidence):
        prob = np.zeros((len(self.states), len(evidence)))
        viterbi = np.zeros((len(self.states), len(evidence)), dtype=int)

        # initializing first column of dynamic array
        for i in range(0, len(self.states)):
            prob[i][0] = np.log(self.emission[x[0]][self.states[i]])

        for i in range(1, len(evidence)):
            for j in range(0, len(self.states)):
                best_prob = 0
                best_state = None
                for k in range(0, len(self.states)):
                    tmp = prob[k][i - 1] * self.transition[self.states[j]][self.states[k]]
                    if tmp > best_prob:
                        best_prob = tmp
                        best_state = k
                prob[j][i] = best_prob * self.emission[evidence[i]][self.states[j]]
                viterbi[j][i] = best_state

        maximum_last_state = np.argmax(prob[:, len(x) - 1])
        output = ''
        for i in range(len(x) - 1, -1, -1):
            output = self.states[maximum_last_state] + output
            maximum_last_state = viterbi[maximum_last_state][i]
        return output


file = open('../Motif Finding/input.txt', 'r')
Lines = file.readlines()
x = Lines[0].strip()
alphabet = Lines[2].split()
states = Lines[4].split()
transition = pd.read_csv(StringIO(''.join(Lines[6:7 + len(states)])), delim_whitespace=True, engine='python')
emission = pd.read_csv(StringIO(''.join(Lines[8 + len(states):9 + len(states) * 2])), delim_whitespace=True,
                       engine='python')
hmm = HMM(emission, transition, alphabet, states)
print(hmm.viterbi(x))

