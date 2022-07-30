#!/usr/bin/env python
import numpy as np
import bioinfo
import gzip
import matplotlib.pyplot as plt
import argparse

def get_args():
    """
    Adds global variables to run different specifications via command line.
    """
    parser = argparse.ArgumentParser(description="Specify parameters")
    parser.add_argument('-r', '--readlen', help='specify read length', type=int)
    parser.add_argument('-f', '--filename', help='specify filename')
    parser.add_argument('-g', '--outputgraph', help='specify output graph name')
    return parser.parse_args()

args = get_args()

def populate_list(file: str) -> tuple[list, int]:
    """Opens a FASTQ file and decodes Phred quality scores to numbers
    accounting for Phred+33. Sums the quality scores for each position
    and counts the total number of lines in the FASTQ file.
    Returns the array and line count."""
    lst = [0] * args.readlen
    with gzip.open(file, 'r') as fq:
        line_count = 0
        for line in fq:
            line = line.decode('ASCII')
            line = line.strip('\n')
            line_count += 1
            # obtain lines with quality scores
            if line_count % 4 == 0:
                # specify for position when converting phred score
                for count, letter in enumerate(line):
                    lst[count] += bioinfo.convert_phred(letter)
    return (lst, line_count)


qscore_list, line_num = populate_list(args.filename)

# calculate mean by dividing the running sum by number of qscore lines
for indice, running_sum in enumerate(qscore_list):
    qscore_list[indice] = running_sum / (line_num / 4)

# create indexes list for easy graphing
indexes = []
for i in enumerate(qscore_list): 
    indexes.append(i[0])

fig, ax = plt.subplots(1, figsize=(20,10))
ax.bar(indexes, qscore_list, width=0.9, color='#87CEEB')    
plt.ylabel("Mean Quality Score", size=15)
plt.xlabel("Base Pair Indice", size=15)
plt.title("Mean Quality Scores for Base Pairs Indices", size=20) 
plt.savefig(args.outputgraph)