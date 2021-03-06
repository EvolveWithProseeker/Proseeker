# Proseeker

# Version:           1.6_LOGOS
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
import shutil
import math
import statistics
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
from scipy.stats import shapiro
from itertools import product
from scipy.stats import kstest
import logomaker as lm
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# DEFINITIONS

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

# DATALOAD

user_input = input("Enter the path to your Proseeker working directory (i.e C:\Proseeker): ")

# user_input = str(sys.argv[1])

collist = ['g', 'k', 'MEG', 'block', 'g1', 'pchoice', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
           'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'truncres', 'cores', 'sites', 'bres', 'fmode', 'fthresh', 'memorylimit',
           'memorythresh','bresnum']
pd.set_option("display.max_colwidth", 10000)
jobstart = pd.read_csv(os.path.join(user_input, 'jobstart.csv'), usecols=collist)
resdir = os.path.join(user_input, datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
bresstr = str(jobstart['bres'])
bresstr = bresstr[5:len(bresstr) - 26]
bresdir = os.path.join(user_input, bresstr)
bresnum = int(jobstart['bresnum'])
generations = int(jobstart['g'])
wide = int(jobstart['k'])
MEG = int(jobstart['MEG'])
block = str(jobstart['block'])
g1 = str(jobstart['g1'])
pchoice = int(jobstart['pchoice'])
trunc = int(jobstart['truncres'])
cores = int(jobstart['cores'])
sites = int(jobstart['sites'])
asslib = str(jobstart["fmode"])
asslibthresh = int(jobstart['fthresh'])
if str(asslib[5:8]) == "YES":
    wide = asslibthresh
memlots = str(jobstart["memorylimit"])
memthresh = int(jobstart['memorythresh'])
d = {}

if pchoice == 1:
    aaprobset = jobstart[
        ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']]
    d['aaprobset'] = aaprobset.values.tolist()
    d['aaprobset'] = [item for sublist in d['aaprobset'] for item in sublist]
else:
    d['aaprobset'] = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                      0.05, 0.05, 0.05, 0.05]

# DIRECTORY_SETUP

os.mkdir(resdir)
shutil.copyfile(os.path.join(user_input, 'jobstart.csv'), os.path.join(resdir, 'used_jobstart.csv'))

#   BRESIMPORT

d1 = {}
d2 = {}

for bs in range(1, sites + 1):
    for x in range(1, bresnum+1):
        bres = pd.read_csv(os.path.join(bresdir, 'p{}.bres{}.csv'.format(x, bs)), header=None, sep=',')
        d['b{}.bres{}.csv'.format(x, bs)] = bres
        d['b{}.bres{}.ind'.format(x, bs)] = list(bres.iloc[:, 13])

        for v in range(0, 13):
            d1['col{}'.format(v)] = list(bres.iloc[:, v])

        for y in range(0, 50):
            set1 = [index for index, element in enumerate(d['b{}.bres{}.ind'.format(x, bs)]) if element == y]
            meanset1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            slopeset1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

            for z in range(0, 13):
                colsel1 = d1['col{}'.format(z)]
                sub1 = [colsel1[i] for i in set1]
                meanset1[z] = sum(sub1) / len(sub1)
                slopex1 = list.copy(sub1)
                for q in range(1, len(slopex1) + 1):
                    slopex1[q - 1] = q
                sub1.sort()
                m, b = np.polyfit(slopex1, sub1, 1)
                slopeset1[z] = m
            d['b{}.bres{}.c{}'.format(x, bs, y)] = meanset1
            d['b{}.cslopes{}.c{}'.format(x, bs, y)] = slopeset1

#   AA values load and normalisation (if the data is pre-normalised it will be unchanged)

aavals = pd.read_csv(os.path.join(user_input, 'ranking.csv'),
                     usecols=['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W',
                              'Y', 'V'],
                     sep=',')

for i in range(0, 544):
    for j in range(0, 20):
        rowmin = min(aavals.iloc[i])
        rowmax = max(aavals.iloc[i])
        val = aavals.iloc[i, j]
        aavals.replace([aavals.iloc[i, j]], (val - rowmin) / (rowmax - rowmin))

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
asslib = str(asslib[5:8])
memlots = str(memlots[5:8])

for g in range(2, generations + 1):
    klstindx = 0

    if asslib != "YES":
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
            del (var)

    else:

        print('Library construction has commenced'.format(g))
        mutinds = sorted(mutinds)
        aaspwr = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
                  'S', 'T', 'W', 'Y', 'V']
        aaspwrprob = d['aaprobset']
        aaspwrsub = []

        for i in range(0, len(aaspwrprob)):
            if aaspwrprob[i] > 0:
                aaspwrsub.append(aaspwr[i])

        comblistmaster = []
        if memlots != "YES":
            for comb in product(aaspwrsub, repeat=len(mutinds)):
                comblistmaster.append(''.join(comb))
        elif memlots == "YES":
            for q in range(0, memthresh):
                selelist = []
                for r in range(0, len(mutinds)):
                    selelist.append(random.choice(aaspwr))
                comblistmaster.append(''.join(selelist))
        random.shuffle(comblistmaster)
        subject = list.copy(g1)
        subject = subject[5:len(subject) - 24]
        for k in range(0, len(comblistmaster)):

            var = list.copy(subject)
            currentcomb = comblistmaster[k]
            for emeg in range(0, len(mutinds)):
                choice = int(mutinds[emeg])
                var[choice] = currentcomb[emeg]

            d['g{}p{}'.format(g, k)] = var
            d['g{}p{}'.format(g, k)] = [item for sublist in d['g{}p{}'.format(g, k)] for item in sublist]
            del var

    # VARIANT ASSESSMENT

    bestmutset = []
    kset = []

    if asslib != "YES":
        print('Gen {} assessment commenced'.format(g))
    else:
        print("Library Assessment Commenced")
        wide = len(comblistmaster)

    for k in range(0, wide):
        if asslib == "YES":
            print("Now assessing variant {} of {}".format(k + 1, wide))
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
            window = [var[x - 6], var[x - 5], var[x - 4], var[x - 3], var[x - 2], var[x - 1], var[x], var[x + 1],
                      var[x + 2], var[x + 3],
                      var[x + 4], var[x + 5], var[x + 6]]

            for v in range(1, bresnum+1):  # For each bres set

                for mbs in range(1, sites + 1):

                    for y in range(0, 50):  # For each cluster
                        set1 = [index for index, element in enumerate(d['b{}.bres{}.ind'.format(v, mbs)]) if
                                element == y]
                        meanset1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                        slopeset1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

                        # Enumeration of each column

                        for z in range(0, 13):
                            colsel = window[z]
                            sub1 = [colsel[i] for i in set1]
                            meanset1[z] = sum(sub1) / len(sub1)
                            slopex1 = list.copy(sub1)
                            for q in range(1, len(slopex1) + 1):
                                slopex1[q - 1] = q
                            sub1.sort()
                            m, b = np.polyfit(slopex1, sub1, 1)
                            slopeset1[z] = m

                        d4['p{}.means{}.c{}.pos{}'.format(k, mbs, y, x)] = [a - b for a, b in
                                                                            zip(d['b{}.bres{}.c{}'.format(v, mbs, y)],
                                                                                meanset1)]
                        d4['p{}.cslopes{}.c{}.pos{}'.format(k, mbs, y, x)] = [a - b for a, b in zip(
                            d['b{}.cslopes{}.c{}'.format(v, mbs, y)], slopeset1)]
                        corr_matrixa = np.corrcoef(d['b{}.bres{}.c{}'.format(v, mbs, y)], meanset1)
                        corra = corr_matrixa[0, 1]
                        d4['p{}.rsq{}.c{}.pos{}'.format(k, mbs, y, x)] = 1 - (corra ** 2)
                        d4['p{}.rmse{}.c{}.pos{}'.format(k, mbs, y, x)] = mean_squared_error(
                            d['b{}.bres{}.c{}'.format(v, mbs, y)], meanset1,
                            squared=False)
                        compaset1 = d['b{}.bres{}.c{}'.format(v, mbs, y)]
                        tempcomp1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

                        for t in range(0, 13):
                            tempcomp1[t] = (compaset1[t] - meanset1[t]) ** 2

                        d4['p{}.sse{}.c{}.pos{}'.format(k, mbs, y, x)] = sum(tempcomp1)

                # BRES SUM from CLUSTERS

                totallist1 = []
                summdif1 = [0, 0, 0, 0, 0]
                ttl = []
                for smbs in range(1, sites + 1):
                    for clust in range(0, 50):
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

        if asslib == "YES":
            if k >= asslibthresh - 1 or k >= wide - 1:
                normtestset = []
                for current in range(0, k):
                    normtestset.append(d['var{}.minscr'.format(current)])
                if wide - 1 >= asslibthresh:
                    if k <=49:
                        stat, p = shapiro(normtestset)
                    else:
                        stat, p = kstest(normtestset,'norm')
                        print('ks stat={} p = {}'.format(stat, p))

                    ## UNIFORM CHECK
                    
                    testmean = statistics.mean(normtestset)
                    testsd = statistics.stdev(normtestset)
                    difmid = min(normtestset) + ((max(normtestset) - min(normtestset))/2)
                    if (difmid - testsd) <= testmean <= (difmid + testsd):
                        p = 1.1
                else:
                    p = 1

                if p > 0.05:
                    if p == 1.1:
                        print("Uniform distribution after {} variants being assessed".format(k + 1))
                    elif p == 1:
                        print("run finalised after {} variants being assessed".format(k + 1))
                    else:
                        print("Normal distribution (p = {}) after {} variants being assessed".format(p, k + 1))
                    
                    print("Finalising run")
                    keepvar = k + 1
                    klist1 = []
                    klist2 = []
                    klist3 = []
                    for j in range(0, k):
                        klist1.append(d['g{}p{}'.format(g, j)])
                        klist2.append(d['var{}.minscr'.format(j)])
                        klist3.append('>g{}p{} '.format(g, j))
                    klist1, klist2, klist3 = zip(*sorted(zip(klist2, klist1, klist3)))
                    with open(os.path.join(resdir, 'LIBRARYvariantscores.tsv'.format(g)), 'w', newline='') as f:
                        writer = csv.writer(f, delimiter='\t')
                        writer.writerows(zip(klist3, klist2, klist1))
                        
                    # Top 25%
                    upper = math.ceil(len(klist1) * 0.25)
                    trunclist = [0] * upper
                    for uprcnt in range(0, upper):
                        trunclist[uprcnt] = klist2[uprcnt]
                    for u in range(0, len(mutinds)):
                        position = mutinds[u]
                        taas = []
                        tcount = []
                        for v in range(0, upper):
                            subject = trunclist[v]
                            isosub = subject[position]
                            if isosub not in taas:
                                taas.append(isosub)
                                tcount.append(1)
                            elif isosub in taas:
                                isopos = taas.index(isosub)
                                tcount[isopos] = tcount[isopos] + 1
                        tperc = list.copy(tcount)
                        for w in range(0, len(tcount)):
                            tperc[w] = tcount[w] / sum(tcount)
                        rows = zip(taas, tcount, tperc)
                        with open(os.path.join(resdir, 'Position_{}_representation_UPPER.tsv'.format(u)), 'w',
                                  newline='') as f:
                            writer = csv.writer(f, delimiter='\t')
                            writer.writerows(rows)

                    toptaas = list.copy(taas)
                    toptcount = list.copy(tcount)
                    toptperc = list.copy(tperc)

                    # Bottom 25%
                    klist2 = list(klist2)
                    klist2.reverse()
                    trunclist = [0] * upper
                    finexcludes = []

                    for uprcnt in range(0, upper):
                        trunclist[uprcnt] = klist2[uprcnt]
                    for u in range(0, len(mutinds)):
                        position = mutinds[u]
                        taas = []
                        tcount = []
                        for v in range(0, upper):
                            subject = trunclist[v]
                            isosub = subject[position]
                            if isosub not in taas:
                                taas.append(isosub)
                                tcount.append(1)
                            elif isosub in taas:
                                isopos = taas.index(isosub)
                                tcount[isopos] = tcount[isopos] + 1
                        tperc = list.copy(tcount)
                        for w in range(0, len(tcount)):
                            tperc[w] = tcount[w] / sum(tcount)

                        rows = zip(taas, tcount, tperc)
                        with open(os.path.join(resdir, 'Position_{}_representation_LOWER.tsv'.format(u)), 'w',
                                  newline='') as f:
                            writer = csv.writer(f, delimiter='\t')
                            writer.writerows(rows)

                        # Intersection
                        
                        subject1 = pd.read_csv(os.path.join(resdir, 'Position_{}_representation_UPPER.tsv'.format(u)),
                                               delimiter='\t', header=None)
                        subject2 = pd.read_csv(os.path.join(resdir, 'Position_{}_representation_LOWER.tsv'.format(u)), delimiter='\t', header=None)

                        del taas
                        del tcount
                        del tperc
                        del toptaas
                        del toptcount
                        del toptperc

                        taas = list(subject2[0])
                        tcount = list(subject2[1])
                        tperc = list(subject2[2])
                        toptaas = list(subject1[0])
                        toptcount = list(subject1[1])
                        toptperc = list(subject1[2])

                        # At this point lists are loaded correctly.

                        for ex in range(0, len(taas)):
                            if taas[ex] in toptaas:
                                ind = toptaas.index(taas[ex])
                                stddev = statistics.stdev(toptcount)
                                if tcount[ex] >= toptcount[ind] and tcount[ex] >= (toptcount[ind] + stddev):
                                    finexcludes.append(
                                        "{} has {}% representation at position {} in the lowest 25% of scores and {}% in "
                                        "the highest 25% of scores and it should be excluded (based on sample size {} of "
                                        "a possible {} combinations).".format(
                                            taas[ex], tperc[ex] * 100, u, toptperc[ind] * 100, keepvar,
                                            len(comblistmaster)))
                                elif tcount[ex] >= toptcount[ind] and tcount[ex] < (toptcount[ind] + stddev):
                                    finexcludes.append(
                                        "{} has {}% representation at position {} in the lowest 25% of scores and {}% in "
                                        "the highest 25% of scores and it could be excluded (based on sample size {} of "
                                        "a possible {} combinations).".format(
                                            taas[ex], tperc[ex] * 100, u, toptperc[ind] * 100, keepvar,
                                            len(comblistmaster)))
                            elif taas[ex] not in toptaas:
                                finexcludes.append(
                                    "{} has {}% representation at position {} in the lowest 25% of scores and 0% in the "
                                    "highest 25% of scores and should be excluded (based on sample size {} of a possible "
                                    "{} combinations).".format(
                                        taas[ex], tperc[ex] * 100, u, keepvar, len(comblistmaster)))

                    dumpfile = open(os.path.join(resdir, "exclusions.txt"), "w")
                    for element in finexcludes:
                        dumpfile.write(element + "\n")
                    dumpfile.close

                    # LOGO SECTION
                    color_scheme = {
                        'A': '#2f4f4f',
                        'R': '#2e8b57',
                        'N': '#800000',
                        'D': '#191970',
                        'C': '#808000',
                        'Q': '#ff0000',
                        'E': '#ff8c00',
                        'G': '#ffd700',
                        'H': '#0000cd',
                        'I': '#ba55d3',
                        'L': '#00ff7f',
                        'K': '#adff2f',
                        'F': '#ffdead',
                        'M': '#ff00ff',
                        'P': '#1e90ff',
                        'S': '#fa8072',
                        'T': '#dda0dd',
                        'W': '#add8e6',
                        'Y': '#ff1493',
                        'V': '#7fffd4'
                    }

                    aalogolis = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

                    for ds in range(0, 3):
                        aaalist = [0] * len(mutinds)
                        aarlist = [0] * len(mutinds)
                        aanlist = [0] * len(mutinds)
                        aadlist = [0] * len(mutinds)
                        aaclist = [0] * len(mutinds)
                        aaqlist = [0] * len(mutinds)
                        aaelist = [0] * len(mutinds)
                        aaglist = [0] * len(mutinds)
                        aahlist = [0] * len(mutinds)
                        aailist = [0] * len(mutinds)
                        aallist = [0] * len(mutinds)
                        aaklist = [0] * len(mutinds)
                        aamlist = [0] * len(mutinds)
                        aaflist = [0] * len(mutinds)
                        aaplist = [0] * len(mutinds)
                        aaslist = [0] * len(mutinds)
                        aatlist = [0] * len(mutinds)
                        aawlist = [0] * len(mutinds)
                        aaylist = [0] * len(mutinds)
                        aavlist = [0] * len(mutinds)

                        aaauprlist = [0] * len(mutinds)
                        aaruprlist = [0] * len(mutinds)
                        aanuprlist = [0] * len(mutinds)
                        aaduprlist = [0] * len(mutinds)
                        aacuprlist = [0] * len(mutinds)
                        aaquprlist = [0] * len(mutinds)
                        aaeuprlist = [0] * len(mutinds)
                        aaguprlist = [0] * len(mutinds)
                        aahuprlist = [0] * len(mutinds)
                        aaiuprlist = [0] * len(mutinds)
                        aaluprlist = [0] * len(mutinds)
                        aakuprlist = [0] * len(mutinds)
                        aamuprlist = [0] * len(mutinds)
                        aafuprlist = [0] * len(mutinds)
                        aapuprlist = [0] * len(mutinds)
                        aasuprlist = [0] * len(mutinds)
                        aatuprlist = [0] * len(mutinds)
                        aawuprlist = [0] * len(mutinds)
                        aayuprlist = [0] * len(mutinds)
                        aavuprlist = [0] * len(mutinds)

                        aaalwrlist = [0] * len(mutinds)
                        aarlwrlist = [0] * len(mutinds)
                        aanlwrlist = [0] * len(mutinds)
                        aadlwrlist = [0] * len(mutinds)
                        aaclwrlist = [0] * len(mutinds)
                        aaqlwrlist = [0] * len(mutinds)
                        aaelwrlist = [0] * len(mutinds)
                        aaglwrlist = [0] * len(mutinds)
                        aahlwrlist = [0] * len(mutinds)
                        aailwrlist = [0] * len(mutinds)
                        aallwrlist = [0] * len(mutinds)
                        aaklwrlist = [0] * len(mutinds)
                        aamlwrlist = [0] * len(mutinds)
                        aaflwrlist = [0] * len(mutinds)
                        aaplwrlist = [0] * len(mutinds)
                        aaslwrlist = [0] * len(mutinds)
                        aatlwrlist = [0] * len(mutinds)
                        aawlwrlist = [0] * len(mutinds)
                        aaylwrlist = [0] * len(mutinds)
                        aavlwrlist = [0] * len(mutinds)

                        for u in range(0, len(mutinds)):
                            subject1 = pd.read_csv(os.path.join(resdir, 'Position_{}_representation_UPPER.tsv'.format(u)),
                                                   delimiter='\t', header=None)
                            subject2 = pd.read_csv(os.path.join(resdir, 'Position_{}_representation_LOWER.tsv'.format(u)), delimiter='\t', header=None)

                            l1 = [0] * len(aalogolis)
                            l2 = [0] * len(aalogolis)
                            l3 = [0] * len(aalogolis)

                            for t in range(0,2):
                                if t == 0:
                                    c1 = list(subject1[0])
                                    c2 = list(subject1[2])
                                elif t == 1:
                                    c1 = list(subject2[0])
                                    c2 = list(subject2[2])

                                for v in range(0, len(aalogolis)):
                                    if aalogolis[v] in c1:
                                        lind = c1.index(aalogolis[v])
                                        if t == 0:
                                            l1[v] = c2[lind]
                                        elif t == 1:
                                            l2[v] = c2[lind]
                            for t in range(0, len(l1)):
                                l3[t] = float(l1[t]) - float(l2[t])

                            aaalist[u] = l3[aalogolis.index('A')]
                            aarlist[u] = l3[aalogolis.index('R')]
                            aanlist[u] = l3[aalogolis.index('N')]
                            aadlist[u] = l3[aalogolis.index('D')]
                            aaclist[u] = l3[aalogolis.index('C')]
                            aaqlist[u] = l3[aalogolis.index('Q')]
                            aaelist[u] = l3[aalogolis.index('E')]
                            aaglist[u] = l3[aalogolis.index('G')]
                            aahlist[u] = l3[aalogolis.index('H')]
                            aailist[u] = l3[aalogolis.index('I')]
                            aallist[u] = l3[aalogolis.index('L')]
                            aaklist[u] = l3[aalogolis.index('K')]
                            aamlist[u] = l3[aalogolis.index('M')]
                            aaflist[u] = l3[aalogolis.index('F')]
                            aaplist[u] = l3[aalogolis.index('P')]
                            aaslist[u] = l3[aalogolis.index('S')]
                            aatlist[u] = l3[aalogolis.index('T')]
                            aawlist[u] = l3[aalogolis.index('W')]
                            aaylist[u] = l3[aalogolis.index('Y')]
                            aavlist[u] = l3[aalogolis.index('V')]

                            aaauprlist[u] = l1[aalogolis.index('A')]
                            aaruprlist[u] = l1[aalogolis.index('R')]
                            aanuprlist[u] = l1[aalogolis.index('N')]
                            aaduprlist[u] = l1[aalogolis.index('D')]
                            aacuprlist[u] = l1[aalogolis.index('C')]
                            aaquprlist[u] = l1[aalogolis.index('Q')]
                            aaeuprlist[u] = l1[aalogolis.index('E')]
                            aaguprlist[u] = l1[aalogolis.index('G')]
                            aahuprlist[u] = l1[aalogolis.index('H')]
                            aaiuprlist[u] = l1[aalogolis.index('I')]
                            aaluprlist[u] = l1[aalogolis.index('L')]
                            aakuprlist[u] = l1[aalogolis.index('K')]
                            aamuprlist[u] = l1[aalogolis.index('M')]
                            aafuprlist[u] = l1[aalogolis.index('F')]
                            aapuprlist[u] = l1[aalogolis.index('P')]
                            aasuprlist[u] = l1[aalogolis.index('S')]
                            aatuprlist[u] = l1[aalogolis.index('T')]
                            aawuprlist[u] = l1[aalogolis.index('W')]
                            aayuprlist[u] = l1[aalogolis.index('Y')]
                            aavuprlist[u] = l1[aalogolis.index('V')]

                            aaalwrlist[u] = l2[aalogolis.index('A')]
                            aarlwrlist[u] = l2[aalogolis.index('R')]
                            aanlwrlist[u] = l2[aalogolis.index('N')]
                            aadlwrlist[u] = l2[aalogolis.index('D')]
                            aaclwrlist[u] = l2[aalogolis.index('C')]
                            aaqlwrlist[u] = l2[aalogolis.index('Q')]
                            aaelwrlist[u] = l2[aalogolis.index('E')]
                            aaglwrlist[u] = l2[aalogolis.index('G')]
                            aahlwrlist[u] = l2[aalogolis.index('H')]
                            aailwrlist[u] = l2[aalogolis.index('I')]
                            aallwrlist[u] = l2[aalogolis.index('L')]
                            aaklwrlist[u] = l2[aalogolis.index('K')]
                            aamlwrlist[u] = l2[aalogolis.index('M')]
                            aaflwrlist[u] = l2[aalogolis.index('F')]
                            aaplwrlist[u] = l2[aalogolis.index('P')]
                            aaslwrlist[u] = l2[aalogolis.index('S')]
                            aatlwrlist[u] = l2[aalogolis.index('T')]
                            aawlwrlist[u] = l2[aalogolis.index('W')]
                            aaylwrlist[u] = l2[aalogolis.index('Y')]
                            aavlwrlist[u] = l2[aalogolis.index('V')]

                        df = pd.DataFrame(list(zip(aaalist, aarlist, aanlist, aadlist, aaclist, aaqlist, aaelist, aaglist, aahlist, aailist, aallist, aaklist, aamlist, aaflist, aaplist, aaslist, aatlist, aawlist, aaylist, aavlist)),
                                          columns=['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
                        plt.ion()
                        logo = lm.Logo(df,
                                           shade_below=.5,
                                           fade_below=.5,
                                           font_name='Courier New',
                                       color_scheme=color_scheme,)

                        logo.style_spines(visible=False)
                        logo.style_spines(spines=['left', 'bottom'], visible=True)
                        logo.style_xticks(rotation=90, fmt='%d', anchor=0)
                        logo.ax.set_ylabel("Difference (Q1-Q4)", labelpad=-1)
                        plt.savefig(os.path.join(resdir, "DifferenceLogo.svg"))
                        plt.close()

                        df = pd.DataFrame(list(
                            zip(aaauprlist, aaruprlist, aanuprlist, aaduprlist, aacuprlist, aaquprlist, aaeuprlist, aaguprlist, aahuprlist,
                                aaiuprlist, aaluprlist, aakuprlist, aamuprlist, aafuprlist, aapuprlist, aasuprlist, aatuprlist, aawuprlist,
                                aayuprlist, aavuprlist)),
                                          columns=['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
                                                   'P', 'S', 'T', 'W', 'Y', 'V'])
                        plt.ion()
                        logo = lm.Logo(df,
                                       shade_below=.5,
                                       fade_below=.5,
                                       font_name='Courier New',
                                       color_scheme=color_scheme)

                        logo.style_spines(visible=False)
                        logo.style_spines(spines=['left', 'bottom'], visible=True)
                        logo.style_xticks(rotation=90, fmt='%d', anchor=0)
                        logo.ax.set_ylabel("Q1 Rep", labelpad=-1)
                        plt.savefig(os.path.join(resdir, "Q1Logo.svg"))
                        plt.close()

                        df = pd.DataFrame(list(zip(aaalwrlist, aarlwrlist, aanlwrlist, aadlwrlist, aaclwrlist, aaqlwrlist, aaelwrlist,
                                                   aaglwrlist, aahlwrlist, aailwrlist, aallwrlist, aaklwrlist, aamlwrlist, aaflwrlist,
                                                   aaplwrlist, aaslwrlist, aatlwrlist, aawlwrlist, aaylwrlist, aavlwrlist)),
                                          columns=['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
                        plt.ion()
                        logo = lm.Logo(df,
                                           shade_below=.5,
                                           fade_below=.5,
                                           font_name='Courier New',
                                       color_scheme=color_scheme,)

                        logo.style_spines(visible=False)
                        logo.style_spines(spines=['left', 'bottom'], visible=True)
                        logo.style_xticks(rotation=90, fmt='%d', anchor=0)
                        logo.ax.set_ylabel("Q4 Rep", labelpad=-1)
                        plt.savefig(os.path.join(resdir, "Q4Logo.svg"))
                        plt.close()
                        exit()

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
        tpctfrc = int(round(wide / trunc))
        klist1 = klist1[0:tpctfrc]
        klist2 = klist2[0:tpctfrc]
        klist2 = map(''.join, klist2)
        klist3 = klist3[0:tpctfrc]

    with open(os.path.join(resdir, 'g{}variantscores.tsv'.format(g)), 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(zip(klist3, klist2, klist1))
