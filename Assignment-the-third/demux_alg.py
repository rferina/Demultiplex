#!/usr/bin/env python
import bioinfo
import argparse
import gzip
import numpy as np
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="Specify parameters")
    parser.add_argument('-r1', '--read1filename', help='specify read 1 filename')
    parser.add_argument('-i1', '--index1filename', help='specify index 1 filename')
    parser.add_argument('-i2', '--index2filename', help='specify index 2 filename')
    parser.add_argument('-r2', '--read2filename', help='specify read 2 filename')
    parser.add_argument('-index', '--validindexfilename', help='specify valid indexes filename')
    parser.add_argument('-outdir', '--outputdirectory', help='specify output directory')
    return parser.parse_args()


args = get_args()

# open index file, place indexes in set
with open(args.validindexfilename) as real_indexes:
    valid_indexes = set()
    real_indexes.readline()
    for line in real_indexes:
        line = line.strip('\n')
        line = str(line)
        columns = line.split()
        valid_indexes.add(columns[4])


def modify_header(line_1, index_1, index_2):
    '''
    Takes in the first line of the record (the header line)
    and index 1 and the reverse complement of index 2.
    Adds index 1 and the reverse compliment of 
    index 2 to the end of the line_1 header. Returns the new header.
    '''
    line_1 = str(line_1)
    labeled_header = line_1 + ' ' + index_1 + '-' + index_2 + '\n'
    return labeled_header


def write_out(filename, line_1, line_2, line_3, line_4):
    '''
    Takes in a filename to write out to, and the
    4 lines of the record that will be written out. 
    Writes out the record to the specified output file.
    Does not return anything.'''
    filename.write(line_1 + line_2 + '\n' + line_3 + '\n' + line_4 + '\n')


# open index-specific output files
files_dict = {}
for index in valid_indexes:
    files_dict[index] = [open(f'{args.outputdirectory}/{index}_r1_matched.fq', 'w'), open(f'{args.outputdirectory}/{index}_r2_matched.fq', 'w')]
    files_dict['unmatched'] = [open(f'{args.outputdirectory}/r1_unmatched.fq', 'w'), open(f'{args.outputdirectory}/r2_unmatched.fq', 'w')]
    files_dict['unknown'] = [open(f'{args.outputdirectory}/r1_unknown.fq', 'w'), open(f'{args.outputdirectory}/r2_unknown.fq', 'w')]


# count the frequencies of matched indexes
index_count_dict = {}

# count the frequencies of unmatched indexes
unmatched_index_dict = {}


# open input FASTQ files, using rt to read the gzipped files
with gzip.open(args.read1filename, 'rt') as read1:
    with gzip.open(args.index1filename, 'rt') as index1:
        with gzip.open(args.index2filename, 'rt') as index2:
            with gzip.open(args.read2filename, 'rt') as read2:
                # initialize counters for unmatched and unknown counts
                unknown_count = 0
                unmatched_count = 0
                #  if the file still has reads, continue looping
                while True:
                    # get a record for read1, index1, index2 and read2 files
                    # get header read1 line, strip newline
                    header_r1 = read1.readline().strip('\n')
                    # get header for index1 file
                    header_i1 = index1.readline().strip('\n')
                    # get header for index2 file
                    header_i2 = index2.readline().strip('\n')
                    # get header read2 line, strip newline
                    header_r2 = read2.readline().strip('\n')
                    # if the header is empty, its the end of the file, so break the while loop
                    if header_r1 == '':
                        break 
                    # if there is a header line, continue reading the record
                    elif header_r1.startswith('@'):
                        # populate arrays with a record, position 0 is header, position 1 is sequence, position 2 is spacer, position 3 is quality scores
                        read1_record = np.array([header_r1.strip('\n'), read1.readline().strip('\n'),  read1.readline().strip('\n'),  read1.readline().strip('\n')])
                        index1_record = np.array([header_i1.strip('\n'), index1.readline().strip('\n'),  index1.readline().strip('\n'),  index1.readline().strip('\n')])
                        index2_record = np.array([header_i2.strip('\n'), index2.readline().strip('\n'),  index2.readline().strip('\n'),  index2.readline().strip('\n')])
                        read2_record = np.array([header_r2.strip('\n'), read2.readline().strip('\n'),  read2.readline().strip('\n'),  read2.readline().strip('\n')])
                        # calculate the reverse complement of index 2
                        rev_comp_index2 = bioinfo.reverse_complement(index2_record[1])
                        # if index 1 is invalid, the read is invalid, write to unknown file
                        if index1_record[1] not in valid_indexes: 
                            r1_indexes_header = modify_header(read1_record[0], index1_record[1], rev_comp_index2)
                            write_out(files_dict['unknown'][0], r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                            r2_indexes_header = modify_header(read2_record[0], index1_record[1], rev_comp_index2)
                            write_out(files_dict['unknown'][1], r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
                            unknown_count += 1
                        # if reverse complement of index 2 is invalid, the read is invalid, write to unknown file
                        elif rev_comp_index2 not in valid_indexes:
                            r1_indexes_header = modify_header(read1_record[0], index1_record[1], rev_comp_index2)
                            write_out(files_dict['unknown'][0], r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                            r2_indexes_header = modify_header(read2_record[0], index1_record[1], rev_comp_index2)
                            write_out(files_dict['unknown'][1], r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
                            unknown_count += 1
                        # if both indexes are valid, see if the indexes pass the quality score cutoff of 30
                        elif index1_record[1] in valid_indexes and rev_comp_index2 in valid_indexes:   
                            # check if index1 doesn't meet quality score cutoff of 30
                            if bioinfo.qual_score(index1_record[3]) < 30:
                                r1_indexes_header = modify_header(read1_record[0], index1_record[1], rev_comp_index2)
                                write_out(files_dict['unknown'][0], r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                                r2_indexes_header = modify_header(read2_record[0], index1_record[1], rev_comp_index2)
                                write_out(files_dict['unknown'][1], r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
                                unknown_count += 1
                            # check if index2 doesn't meet quality score cutoff of 30
                            elif bioinfo.qual_score(index2_record[3]) < 30:
                                r1_indexes_header = modify_header(read1_record[0], index1_record[1], rev_comp_index2)
                                write_out(files_dict['unknown'][0], r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                                r2_indexes_header = modify_header(read2_record[0], index1_record[1], rev_comp_index2)
                                write_out(files_dict['unknown'][1], r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
                                unknown_count += 1
                            # now that both indexes are equal to or above quality score cutoff, see if index 1 matches rev comp index 2               
                            # if not matching, write to unmatched
                            elif index1_record[1] != rev_comp_index2:
                                r1_indexes_header = modify_header(read1_record[0], index1_record[1], rev_comp_index2)
                                write_out(files_dict['unmatched'][0], r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])                                    
                                r2_indexes_header = modify_header(read2_record[0], index1_record[1], rev_comp_index2)
                                write_out(files_dict['unmatched'][1], r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
                                unmatched_count += 1
                                # count how many mismatches
                                if index1_record[1] in unmatched_index_dict:
                                    unmatched_index_dict[index1_record[1]] += 1
                                else:
                                    unmatched_index_dict[index1_record[1]] = 1
                            # indexes are matching; passed quality score cutoff
                            elif index1_record[1] == rev_comp_index2:           
                                for a_key in files_dict:
                                    if index1_record[1] == a_key:
                                        r1_indexes_header = modify_header(read1_record[0], index1_record[1], rev_comp_index2)
                                        write_out(files_dict[index1_record[1]][0], r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                                        r2_indexes_header = modify_header(read2_record[0], index1_record[1], rev_comp_index2)
                                        write_out(files_dict[index1_record[1]][1], r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
                                # count the frequency of matching indexes
                                if index1_record[1] in index_count_dict:
                                    index_count_dict[index1_record[1]] += 1
                                else:
                                    index_count_dict[index1_record[1]] = 1
                        
                            
# find the total number of matched indexes 
total_matched = sum(index_count_dict.values())
# find the total reads in the file
total_count = total_matched + unmatched_count + unknown_count
index_percentages = {}
# calculate percentages of matched indexes over all reads
for index in index_count_dict: 
    index_percentages[index] = index_count_dict[index] / total_count * 100
# Return unmatched and unknown counts, matched dictionary counts, stats                           
print('Number of unmatched records:', unmatched_count)
print('Number of unknown records:', unknown_count) 
print('Number of total matched records:', total_matched)
print('Total Number of Records:', total_count)
print('Unmatched index counts:', unmatched_index_dict)
print('Matched index counts:', index_count_dict)
print('Matched index percentages:', index_percentages)


# graph pie chart with percentages of reads from each matched index
index_perc = list(index_percentages.values())
index_labels = list(index_percentages.keys())
plt.pie(index_perc, labels=index_labels, autopct='%1.1f%%', rotatelabels=True)
plt.title("Percentages of Valid Reads from Each Index", size=15, pad=28) 
plt.savefig('valid_reads.png')

# graph frequency distribution
read_type = ['Matched', 'Unmatched', 'Unknown']
read_counts = [total_matched, unmatched_count, unknown_count]
fig, ax = plt.subplots(1, figsize=(15,10))
ax.bar(read_type, read_counts, width=0.9, color='#87CEEB')    
plt.ylabel("Frequency", size=18)
plt.yticks(fontsize=15)
plt.xlabel("Read Type", size=18)
plt.xticks(fontsize=15)
plt.title("Frequency of Demultiplexed Read Types", size=20) 
plt.savefig('frequency.png')

# graph index hopping distribution
unmatch_index_keys = list(unmatched_index_dict.keys())
unmatch_values = list(unmatched_index_dict.values())
fig, ax = plt.subplots(1, figsize=(15,10))
ax.bar(unmatch_index_keys, unmatch_values, width=0.9, color='#87CEEB')    
plt.ylabel("Frequency", size=18)
plt.yticks(fontsize=15)
plt.xlabel("Index", size=18)
plt.xticks(fontsize=15, rotation=45, ha="right")
plt.title("Frequency of Index Hopping", size=20) 
plt.savefig('hopping_dist.png')

# graph index matching distribution
match_index_keys = list(index_count_dict.keys())
match_values = list(index_count_dict.values())
fig, ax = plt.subplots(1, figsize=(15,10))
ax.bar(match_index_keys, match_values, width=0.9, color='#87CEEB')    
plt.ylabel("Frequency", size=18)
plt.yticks(fontsize=15)
plt.xlabel("Index", size=18)
plt.xticks(fontsize=15, rotation=45, ha="right")
plt.title("Frequency of Matched Indexes", size=20) 
plt.savefig('matching_dist.png')

#for loop to close matched output files             
for filename in files_dict.values():
    filename[0].close()
    filename[1].close()




