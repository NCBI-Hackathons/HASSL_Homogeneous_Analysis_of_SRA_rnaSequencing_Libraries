from collections import defaultdict, Counter
import csv
import itertools

import numpy as np
import matplotlib as mpl
mpl.use('pdf')  # enables using matplotpib over an ssh connection without X
import matplotlib.pyplot as plt
from sqt.io.fasta import FastaReader, FastaWriter, FastqReader,\
    FastqWriter, SequenceReader


# we can use standard Python functions and call them from different rules below
def readlength_histogram(filename, offset=0):
    """return a Counter representing the histogram of read length in file <filename>"""
    hist = Counter()
    with SequenceReader(filename) as fr:
        for record in fr:
            hist[len(record.sequence) + offset] += 1
    return hist

def write_counter_to_csv(counter, filename, delimiter="\t"):
    with open(filename, "wt") as f:
        out = csv.writer(f, delimiter=delimiter, lineterminator="\n")
        for l in range(max(counter)+1): 
            out.writerow((l, counter[l]))
        del out

def plot_histogram(hist, outname, xlabel="", ylabel="frequency"):
    """use matplotlib.pyplot to plot an m*2-matrix as histogram and save it to a file"""
    plt.bar(hist[:,0], hist[:,1])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(outname)
    plt.close()

def plot_counter(counter, outname, xlabel="", ylabel="frequency"):
    """use matplotlib.pyplot to plot a Counter as histogram and save it to a file"""
    x = sorted(counter.keys())
    y = [counter[k] for k in x]
    plt.bar(x,y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(outname)
    plt.close()

