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
    # parser.add_argument('-o', '--outputfilename', help='specify output filename')
    parser.add_argument('-g', '--outputgraph', help='specify output graph name')
    # parser.add_argument('-r', '--resultsfilename', help='specify output results filename')
    return parser.parse_args()


args = get_args()

# PS4

def populate_list(file: str) -> tuple[list, int]:
    """Opens a FASTQ file and decodes Phred quality scores to numbers
    accounting for Phred+33. Sums the quality scores for each position
    and counts the total number of lines in the FASTQ file.
    Returns the array and line count."""

    # lst = bioinfo.init_list(args.readlen)
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

for indice, running_sum in enumerate(qscore_list):
    qscore_list[indice] = running_sum / (line_num / 4)
    output_list = str(indice)  +'\t' + str(qscore_list[indice])
    print(output_list)

# print(qscore_list)
# read1_list, read1_line_num = populate_list('/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz')
# read1_list, read1_line_num = populate_list('/projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-first/test.fq.gz')
# print(read1_list)

# # index1_list, index1_line_num = bioinfo.populate_list()
# # index2_list, index2_line_num = bioinfo.populate_list()
# # read2_list, read2_line_num = bioinfo.populate_list()




# PS9

# def populate_array(file): 
#     """Opens a FASTQ file and decodes Phred quality scores to numbers
#     accounting for Phred+33. Sums the quality scores for each position
#     and counts the total number of lines in the FASTQ file.
#     Returns the array and line count."""
#     # with gzip.open(file, 'r') as fq:
#     #     # initialize qscores array
#     #     line_count = 0
#     #     for line in fq:
#     #         line_count += 1
#     #     qscores = np.zeros((101, line_count//4), dtype=np.int64)
    
#     # qscores = [0] * args.read_len
#     with gzip.open(file, 'r') as fq:
#         pos = 0
#         line_count = 0
#         for line in fq:
#             line = line.decode('ASCII')
#             line = line.strip('\n')
#             line_count += 1
#             # obtain lines with quality scores
#             if line_count % 4 == 0:
#                 # specify for position when converting phred score
#                 for letter in range(len(line)):
#                     qscores[letter, pos] = bioinfo.convert_phred(line[letter])
#                 pos +=1
#     return (qscores, line_count)

# read1_list, read1_line_num = populate_array('/projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-first/test.fq.gz')

# qscore_list, line_num = populate_array(args.filename)

# print(read1_list)

# mean = np.mean(qscore_list, axis=1)
# mean = np.mean(qscore_list)
# print(mean)
# print(len(mean))


# for index in range(len(mean)):
#     print(str(index) + '\t' + str(mean[index]))


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