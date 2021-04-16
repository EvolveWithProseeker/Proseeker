# Simplified Bres Maker

# Version: 1.0

#Python Version: 2.0

# IMPORTS

import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from numpy import asarray
from numpy import savetxt
import sys
import os

# DEFINITIONS

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

# DATALOAD

#user_input = str(sys.argv[1])
#ranking = str(sys.argv[2])
#working = str(sys.argv[3])
#iterations = int(sys.argv[4])
#trys = int(sys.argv[5])

user_input = "D:/Proseeker/PSEUDOTALEDEETS.csv"
ranking = "D:/Proseeker/ranking.csv"
working = "D:/Proseeker"
iterations = 1000000
trys = 1000

aavals = pd.read_csv(ranking, usecols=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'],
                     sep =',')

d = {}

for i in range(0,544):
    for j in range(0,20):
        rowmin = min(aavals.iloc[i])
        rowmax = max(aavals.iloc[i])
        val = aavals.iloc[i, j]
        aavals.replace([aavals.iloc[i, j]], (val - rowmin)/(rowmax - rowmin))

d['A'] = list(aavals['A'])
d['R'] = list(aavals['R'])
d['D'] = list(aavals['D'])
d['N'] = list(aavals['N'])
d['C'] = list(aavals['C'])
d['E'] = list(aavals['E'])
d['Q'] = list(aavals['Q'])
d['G'] = list(aavals['G'])
d['H'] = list(aavals['H'])
d['I'] = list(aavals['I'])
d['L'] = list(aavals['L'])
d['K'] = list(aavals['K'])
d['M'] = list(aavals['M'])
d['F'] = list(aavals['F'])
d['P'] = list(aavals['P'])
d['S'] = list(aavals['S'])
d['T'] = list(aavals['T'])
d['W'] = list(aavals['W'])
d['Y'] = list(aavals['Y'])
d['V'] = list(aavals['V'])


library = pd.read_csv(user_input, header=None, sep=',')

seqs = library[0]
sites = library[1]

# PROCESSING

for x in range(0, len(seqs)):
    subjectstd = list(seqs[x])
    subject = list.copy(subjectstd)

    for p in range(0,len(subjectstd)):
        subject.append(subjectstd[p])

    for z in range(0, len(subject)):
        if subject[z] == 'A':
            subject[z] = d['A']
        elif subject[z] == 'a':
            subject[z] = d['A']
        elif subject[z] == 'R':
            subject[z] = d['R']
        elif subject[z] == 'r':
            subject[z] = d['R']
        elif subject[z] == 'N':
            subject[z] = d['N']
        elif subject[z] == 'n':
            subject[z] = d['N']
        elif subject[z] == 'D':
            subject[z] = d['D']
        elif subject[z] == 'd':
            subject[z] = d['D']
        elif subject[z] == 'C':
            subject[z] = d['C']
        elif subject[z] == 'c':
            subject[z] = d['C']
        elif subject[z] == 'Q':
            subject[z] = d['Q']
        elif subject[z] == 'q':
            subject[z] = d['Q']
        elif subject[z] == 'E':
            subject[z] = d['E']
        elif subject[z] == 'e':
            subject[z] = d['E']
        elif subject[z] == 'G':
            subject[z] = d['G']
        elif subject[z] == 'g':
            subject[z] = d['G']
        elif subject[z] == 'H':
            subject[z] = d['H']
        elif subject[z] == 'h':
            subject[z] = d['H']
        elif subject[z] == 'I':
            subject[z] = d['I']
        elif subject[z] == 'i':
            subject[z] = d['I']
        elif subject[z] == 'L':
            subject[z] = d['L']
        elif subject[z] == 'l':
            subject[z] = d['L']
        elif subject[z] == 'K':
            subject[z] = d['K']
        elif subject[z] == 'k':
            subject[z] = d['K']
        elif subject[z] == 'M':
            subject[z] = d['M']
        elif subject[z] == 'm':
            subject[z] = d['M']
        elif subject[z] == 'F':
            subject[z] = d['F']
        elif subject[z] == 'f':
            subject[z] = d['F']
        elif subject[z] == 'P':
            subject[z] = d['P']
        elif subject[z] == 'p':
            subject[z] = d['P']
        elif subject[z] == 'S':
            subject[z] = d['S']
        elif subject[z] == 's':
            subject[z] = d['S']
        elif subject[z] == 'T':
            subject[z] = d['T']
        elif subject[z] == 't':
            subject[z] = d['T']
        elif subject[z] == 'W':
            subject[z] = d['W']
        elif subject[z] == 'w':
            subject[z] = d['W']
        elif subject[z] == 'Y':
            subject[z] = d['Y']
        elif subject[z] == 'y':
            subject[z] = d['Y']
        elif subject[z] == 'V':
            subject[z] = d['V']
        elif subject[z] == 'v':
            subject[z] = d['V']

    subjectsites = str(sites[x])
    splits = find(subjectsites, ':')
    splits.append(len(subjectsites))
    if sum(splits) > 0:
        for q in range(len(splits)):
            if q == 0:
                subpos = int(subjectsites[0:splits[q]])
            else:
                subpos = int(subjectsites[splits[q-1]+1:splits[q]])

            breswindow = list((subject[subpos-6], subject[subpos-5], subject[subpos-4], subject[subpos-3],
                                         subject[subpos-2], subject[subpos-1], subject[subpos], subject[subpos+1],
                                         subject[subpos+2], subject[subpos+3], subject[subpos+4], subject[subpos+5],
                                         subject[subpos+6]))

            breswindow = np.column_stack(breswindow)
            kmeans = KMeans(n_clusters=50, n_init=trys, max_iter=iterations, algorithm="full")
            kmeans.fit(breswindow)
            clusters = kmeans.labels_
            breswindow = np.insert(breswindow, 13, clusters, axis=1)
            savetxt(os.path.join(working, 'p{}.bres{}.csv'.format(x+1, q+1)), breswindow, delimiter=',', fmt='%f')
