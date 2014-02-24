#!/usr/bin/env python

###############################################################################
#
# deltagc_repeat - calculate and plot the difference of GC ratio between the
#                  siRNA-accumulated regions and the non-siRNA-accumulated
#                  regions, and compare with the randomized tandem repeats
#
# written by Wen-Wei Liao (wwliao@gate.sinica.edu.tw)
#
###############################################################################

from __future__ import division
import argparse
import cPickle as pickle
from random import shuffle
from collections import defaultdict
from Bio import SeqIO
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

def overlap(a, b):
    """
    Compute the number of overlapping bases between two regions
    """
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def gc_ratio(seq):
    """
    Compute the GC ratio of a given DNA sequence
    """
    gc_count = seq.count('G') + seq.count('C')
    at_count = seq.count('A') + seq.count('T')
    return gc_count / (at_count + gc_count)

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('sirna')
    parser.add_argument('repeat')
    parser.add_argument('ref')
    return parser

def deltagc(sirna, trepeat, ref, rand):
    rna_regions = defaultdict(list)
    with open(sirna) as infile:
        for i, line in enumerate(infile):
            if i != 0:
                line = line.strip().split('\t')
                l = int(line[1])
                chr = line[2]
                strand = line[3]
                pos = int(line[4])
                if strand == 'w':
                    rna_regions[chr].append((pos - 1, pos + l - 1))
                elif strand == 'c':
                    rna_regions[chr].append((pos - l, pos))

    repeat = defaultdict(list)
    with open(trepeat) as infile:
        for i, line in enumerate(infile):
            if i != 0:
                line = line.strip().split('\t')
                chr = line[0]
                start = int(line[1]) - 1    # change to 0-based coordinate
                end = int(line[2])
                repeat[chr].append((start, end))

    rna_repeat = defaultdict(lambda: defaultdict(list))
    for chr in repeat:
        for rna in rna_regions[chr]:
            for rp in repeat[chr]:
                if overlap(rna, rp):
                    rna_repeat[chr][rp].append((max(0, rna[0]-rp[0]), min(rp[1]-rp[0], rna[1]-rp[0])))
                    break

    rna_flag = {}
    for chr in rna_repeat:
        rna_flag[chr] = {}
        for rp in rna_repeat[chr]:
            rna_flag[chr][rp] = [0]*(rp[1]-rp[0])
            for rna in rna_repeat[chr][rp]:
                for i in xrange(rna[0], rna[1]):
                    rna_flag[chr][rp][i] = 1

    with open(ref) as infile:
        ref = SeqIO.to_dict(SeqIO.parse(infile, 'fasta'))
        for chr in ref:
            ref[chr] = str(ref[chr].seq).upper()

    coverall = 0
    deltagc = []
    for chr in rna_flag:
        for rp in rna_flag[chr]:
            wrna = []
            orna = []
            seq = ref[chr][rp[0]:rp[1]]
            if rand:
                seq = list(seq)
                shuffle(seq)
                seq = ''.join(seq)
            for flag, base in zip(rna_flag[chr][rp], seq):
                if flag == 1:
                    wrna.append(base)
                else:
                    orna.append(base)
            if wrna and orna:
                deltagc.append(gc_ratio(wrna) - gc_ratio(orna))
            else:
                coverall += 1
    return deltagc

def boxplot(real, rand):
    colors = {'green': (143/255, 188/255,  87/255),
               'blue': (112/255, 193/255, 222/255),
                'red': (226/255, 121/255,  76/255),
              'black': ( 35/255,  31/255,  30/255),
               'gray': (125/255, 125/255, 125/255),
              'white': (252/255, 244/255, 233/255)}

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_axes([0.15, 0.15, 0.75, 0.75], axisbg=colors['white'])
    data =[real, rand]
    bp = ax.boxplot(data, widths=0.4, patch_artist=True)
    bp['boxes'][0].set_color(colors['blue'])
    bp['boxes'][1].set_color(colors['blue'])
    plt.setp(bp['boxes'], linewidth=1.5)
    plt.setp(bp['caps'], color=colors['black'])
    plt.setp(bp['caps'], linewidth=1.5)
    plt.setp(bp['medians'], color=colors['black'])
    plt.setp(bp['medians'], linewidth=1.5)
    plt.setp(bp['whiskers'], color=colors['black'])
    plt.setp(bp['whiskers'], linestyle='-')
    plt.setp(bp['whiskers'], linewidth=1.5)
    plt.setp(bp['whiskers'], zorder=-1)
    plt.setp(bp['fliers'], color=colors['gray'])
    plt.setp(bp['fliers'], marker='x')
    plt.setp(bp['fliers'], markeredgewidth=1.5)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['left'].set_linewidth(3)
    ax.tick_params(direction='out', width=3, length=5, labelsize='x-large', top='off', bottom='off', right='off')
    ax.set_xticklabels(['Tandem Repeats\nwith siRNA', 'Randomized\nTandem Repeats'])
    for label in ax.get_xticklabels():
        label.set_fontweight('bold')
    for label in ax.get_yticklabels():
        label.set_fontweight('bold')
    ax.set_ylabel('Difference of GC Ratio', size='x-large', weight='bold')
    ax.axhline(linestyle='--', linewidth=1, color='k', zorder=-1)
    plt.savefig('deltagc_boxplot.png', dpi=300)
    plt.savefig('deltagc_boxplot.eps')
    plt.close(fig)

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    real = deltagc(args.sirna, args.repeat, args.ref, False)
    rand = deltagc(args.sirna, args.repeat, args.ref, True)
    t, p = ttest_ind(real, rand, equal_var=False)
    print 't-statistic: {0}'.format(t)
    print '    p-value: {0}'.format(p)
    boxplot(real, rand)
