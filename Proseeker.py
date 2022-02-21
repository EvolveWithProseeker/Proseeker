# Proseeker

# Version:           1.0
# Python Version:    3.9

# IMPORTS

import sys
import subprocess
import pkg_resources
import os
from datetime import datetime
import random
import warnings
import csv

required = {'numpy', 'pandas', 'sklearn'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

if missing:
    python = sys.executable
    subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)

import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error

warnings.filterwarnings("ignore")

# DEFINITIONS

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

# DATALOAD

user_input = input("Enter the path to your Proseeker working directory (i.e C:\Proseeker): ")

# user_input = str(sys.argv[1])

collist = ['g', 'k', 'MEG', 'block', 'g1', 'pchoice', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
           'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'truncres', 'cores', 'sites', 'bres']

jobstart = pd.read_csv(os.path.join(user_input, 'jobstart.csv'), usecols=collist)
resdir = os.path.join(user_input, datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
bresstr = str(jobstart['bres'])
bresstr = bresstr[5:len(bresstr)-26]
bresdir = os.path.join(user_input, bresstr)
generations = int(jobstart['g'])
wide = int(jobstart['k'])
MEG = int(jobstart['MEG'])
block = str(jobstart['block'])
g1 = str(jobstart['g1'])
pchoice = int(jobstart['pchoice'])
trunc = int(jobstart['truncres'])
cores = int(jobstart['cores'])
sites = int(jobstart['sites'])

d = {}

if pchoice == 1:
    aaprobset = jobstart[['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']]
    d['aaprobset'] = aaprobset.values.tolist()
    d['aaprobset'] = [item for sublist in d['aaprobset'] for item in sublist]
else:
    d['aaprobset'] = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                      0.05, 0.05, 0.05, 0.05]

# DIRECTORY_SETUP

os.mkdir(resdir)

#   BRESIMPORT

d1 = {}
d2 = {}

for bs in range(1, sites +1):
    for x in range(1, 31):
        bres = pd.read_csv(os.path.join(bresdir, 'p{}.bres{}.csv'.format(x, bs)), header=None, sep=',')
        d['b{}.bres{}.csv'.format(x, bs)] = bres
        d['b{}.bres{}.ind'.format(x, bs)] = list(bres.iloc[:, 13])

        for v in range(0,13):
            d1['col{}'.format(v)] = list(bres.iloc[:, v])

        for y in range (1, 51):
            set1 = [index for index, element in enumerate(d['b{}.bres{}.ind'.format(x, bs)]) if element == y]
            meanset1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            slopeset1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

            for z in range(0, 13):
                colsel1 = d1['col{}'.format(z)]
                sub1 = [colsel1[i] for i in set1]
                meanset1[z] = sum(sub1) / len(sub1)
                slopex1 = list.copy(sub1)
                for q in range(1, len(slopex1) + 1):
                    slopex1[q-1] = q
                sub1.sort()
                m, b = np.polyfit(slopex1, sub1, 1)
                slopeset1[z] = m
            d['b{}.bres{}.c{}'.format(x, bs, y)] = meanset1
            d['b{}.cslopes{}.c{}'.format(x, bs, y)] = slopeset1

 #   AA values load and normalisation (if the data is pre-normalised it will be unchanged)

aavals = pd.read_csv(os.path.join(user_input, 'ranking.csv'),
                     usecols=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'],
                     sep =',')

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

mutinds = find(block, 'M')
mutinds = [x - 5 for x in mutinds]
readorder = list.copy(mutinds)

# MUTAGENESIS

g1 = list(g1)
bestmutset = []

for g in range(2, generations+1):
    klstindx = 0


    print('Gen {} mutagenesis commenced'.format(g))


    for k in range(0, wide):

        if g == 2:
            subject = list.copy(g1)
            subject = subject[5:len(subject) - 24]
        else:
            if k / 10 >= klstindx + 1:
                klstindx += 1
            subject = bestmutset[klstindx]

        var = list.copy(subject)

        for emeg in range(0, MEG):

            random.shuffle(mutinds)

            mutchoice = np.random.choice(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
                                          'S', 'T', 'W', 'Y', 'V'], 1, p=list(d['aaprobset']))

            choice = int(mutinds[0])
            var[choice] = mutchoice

        d['g{}p{}'.format(g, k)] = var
        d['g{}p{}'.format(g, k)] = [item for sublist in d['g{}p{}'.format(g, k)] for item in sublist]

        del(var)

# VARIANT ASSESSMENT

    bestmutset = []
    kset = []

    print('Gen {} assessment commenced'.format(g))

    for k in range(0, wide):
        var = list(d['g{}p{}'.format(g, k)])

        for res in range(0, len(var)):
            if var[res] == 'A':
                var[res] = d['A']
            elif var[res] == 'a':
                var[res] = d['A']
            elif var[res] == 'R':
                var[res] = d['R']
            elif var[res] == 'r':
                var[res] = d['R']
            elif var[res] == 'N':
                var[res] = d['N']
            elif var[res] == 'n':
                var[res] = d['N']
            elif var[res] == 'D':
                var[res] = d['D']
            elif var[res] == 'd':
                var[res] = d['D']
            elif var[res] == 'C':
                var[res] = d['C']
            elif var[res] == 'c':
                var[res] = d['C']
            elif var[res] == 'Q':
                var[res] = d['Q']
            elif var[res] == 'q':
                var[res] = d['Q']
            elif var[res] == 'E':
                var[res] = d['E']
            elif var[res] == 'e':
                var[res] = d['E']
            elif var[res] == 'G':
                var[res] = d['G']
            elif var[res] == 'g':
                var[res] = d['G']
            elif var[res] == 'H':
                var[res] = d['H']
            elif var[res] == 'h':
                var[res] = d['H']
            elif var[res] == 'I':
                var[res] = d['I']
            elif var[res] == 'i':
                var[res] = d['I']
            elif var[res] == 'L':
                var[res] = d['L']
            elif var[res] == 'l':
                var[res] = d['L']
            elif var[res] == 'K':
                var[res] = d['K']
            elif var[res] == 'k':
                var[res] = d['K']
            elif var[res] == 'M':
                var[res] = d['M']
            elif var[res] == 'm':
                var[res] = d['M']
            elif var[res] == 'F':
                var[res] = d['F']
            elif var[res] == 'f':
                var[res] = d['F']
            elif var[res] == 'P':
                var[res] = d['P']
            elif var[res] == 'p':
                var[res] = d['P']
            elif var[res] == 'S':
                var[res] = d['S']
            elif var[res] == 's':
                var[res] = d['S']
            elif var[res] == 'T':
                var[res] = d['T']
            elif var[res] == 't':
                var[res] = d['T']
            elif var[res] == 'W':
                var[res] = d['W']
            elif var[res] == 'w':
                var[res] = d['W']
            elif var[res] == 'Y':
                var[res] = d['Y']
            elif var[res] == 'y':
                var[res] = d['Y']
            elif var[res] == 'V':
                var[res] = d['V']
            elif var[res] == 'v':
                var[res] = d['V']

        # FORMULATION OF DESCRIPTIVE STATISTICS

        d4 = {}

        vartot = []

        var0 = var[0]
        var1 = var[1]
        var2 = var[2]
        var3 = var[3]
        var4 = var[4]
        var5 = var[5]

        var.append(var0)
        var.append(var1)
        var.append(var2)
        var.append(var3)
        var.append(var4)
        var.append(var5)

        for x in readorder:
            window = [var[x-6], var[x-5], var[x-4], var[x-3], var[x-2], var[x-1], var[x], var[x+1], var[x+2], var[x+3],
                      var[x+4], var[x+5], var[x+6]]

            for v in range(1, 31): # For each bres set

                for mbs in range(1, sites + 1):

                    for y in range(1, 51): # For each cluster
                        set1 = [index for index, element in enumerate(d['b{}.bres{}.ind'.format(v, mbs)]) if element == y]
                        meanset1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                        slopeset1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

                        # Enumeration of each column

                        for z in range(0, 13):
                            colsel = window[z]
                            sub1 = [colsel[i] for i in set1]
                            meanset1[z] = sum(sub1) / len(sub1)
                            slopex1 = list.copy(sub1)
                            for q in range(1, len(slopex1) + 1):
                                slopex1[q-1] = q
                            sub1.sort()
                            m, b = np.polyfit(slopex1, sub1, 1)
                            slopeset1[z] = m

                        d4['p{}.means{}.c{}.pos{}'.format(k, mbs, y, x)] = [a - b for a, b in zip(d['b{}.bres{}.c{}'.format(v, mbs, y)], meanset1)]
                        d4['p{}.cslopes{}.c{}.pos{}'.format(k, mbs, y, x)] = [a - b for a, b in zip(d['b{}.cslopes{}.c{}'.format(v, mbs, y)], slopeset1)]
                        corr_matrixa = np.corrcoef(d['b{}.bres{}.c{}'.format(v, mbs, y)], meanset1)
                        corra = corr_matrixa[0, 1]
                        d4['p{}.rsq{}.c{}.pos{}'.format(k, mbs, y, x)] = 1- (corra ** 2)
                        d4['p{}.rmse{}.c{}.pos{}'.format(k, mbs, y, x)] = mean_squared_error(d['b{}.bres{}.c{}'.format(v, mbs, y)], meanset1,
                                                                                squared=False)
                        compaset1 = d['b{}.bres{}.c{}'.format(v, mbs, y)]
                        tempcomp1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

                        for t in range(0,13):
                            tempcomp1[t] = (compaset1[t] - meanset1[t]) ** 2

                        d4['p{}.sse{}.c{}.pos{}'.format(k, mbs, y, x)] = sum(tempcomp1)

                # BRES SUM from CLUSTERS

                totallist1 = []
                summdif1 = [0, 0, 0, 0, 0]
                ttl = []
                for smbs in range(1, sites + 1):
                    for clust in range(1, 51):

                        summdif1[0] = sum(d4['p{}.means{}.c{}.pos{}'.format(k, smbs, clust, x)])
                        summdif1[1] = sum(d4['p{}.cslopes{}.c{}.pos{}'.format(k, smbs, clust, x)])
                        summdif1[2] = d4['p{}.rsq{}.c{}.pos{}'.format(k, smbs, clust, x)]
                        summdif1[3] = d4['p{}.rmse{}.c{}.pos{}'.format(k, smbs, clust, x)]
                        summdif1[4] = d4['p{}.sse{}.c{}.pos{}'.format(k, smbs, clust, x)]
                        total1 = sum(summdif1)
                        totallist1.append(total1)
                        tsum = sum(totallist1)
                        ttl.append(tsum)

                vartot.append(sum(ttl))

        d['var{}.minscr'.format(k)] = min(vartot)

    klist1 = []
    klist2 = []
    klist3 = []

    for k in range(0, wide):
        klist1.append(d['g{}p{}'.format(g, k)])
        klist2.append(d['var{}.minscr'.format(k)])
        klist3.append('>g{}p{} '.format(g, k))

    klist1, klist2, klist3 = zip(*sorted(zip(klist2, klist1, klist3)))

    bestmutset = klist2

    if trunc > 0:
        tpctfrc = int(round(wide/trunc))
        klist1 = klist1[0:tpctfrc]
        klist2 = klist2[0:tpctfrc]
        klist2 = map(''.join, klist2)
        klist3 = klist3[0:tpctfrc]

    with open(os.path.join(resdir, 'g{}variantscores.tsv'.format(g)), 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(klist3, klist2, klist1))
