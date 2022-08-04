Rachel Ferina 
Lab Notebook: Demux Part 1

Python version: Python 3.10.4
Environment: No environment used

**26 July 2022**

Created public Demultiplex repository by using the template from https://github.com/Leslie-C/Demultiplexing. 
Cloned Demultiplex repository onto talapas in the directory demux that I created, /projects/bgmp/rferina/bioinfo/Bi622/demux

The data is located on talapas in /projects/bgmp/shared/2017_sequencing/
It should not be unzipped or copied, as the files are very large (300,000,000+ lines).
    1294_S1_L008_R1_001.fastq.gz
    1294_S1_L008_R2_001.fastq.gz
    1294_S1_L008_R3_001.fastq.gz
    1294_S1_L008_R4_001.fastq.gz

Worked on Part 1:

``` zcat 1294_S1_L008_R1_001.fastq.gz | head -25```
``` zcat 1294_S1_L008_R2_001.fastq.gz | head -25```
``` zcat 1294_S1_L008_R3_001.fastq.gz | head -25```
``` zcat 1294_S1_L008_R4_001.fastq.gz | head -25```
| File name | Contents | 
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 | 
| 1294_S1_L008_R4_001.fastq.gz | Read 2 |


Read length:
```zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc```
output:
```
 1       1     102
```
The read length for reads is 102 - 1, as the line as a newline character, so read length is 101.


``` zcat 1294_S1_L008_R2_001.fastq.gz | head -2 | tail -1 | wc```
output:
```
 1       1       9
 ```
The read length for the indexes is 9 - 1, so read length is 8.


If the data was Phred+64 encoded, the quality score lines would include lowercase letters indicating high decimal values, as when converting ABSCII to Phred score, the values would be high enough to subtract 64 from. 
I searched for lowercase letters in the quality score lines on the index file to reduce runtime as they have shorter sequences.
```
zcat 1294_S1_L008_R2_001.fastq.gz | sed -n 4~4p | grep -E "[a-z]+"
```

No output was found, so this data is Phred+33 encoded.

I created the python script demux_processing.py to generate histograms of mean quality scores.

When trying to make a distribution of the mean quality scores, I considered using the strategy with a list from PS4 and the strategy using an array from PS9.


I thought about using PS9, but I had hard-coded the read lengths. I tested if numpy would work to calculate the mean on a list, and it did.
    lst = [4, 5, 7, 8, 9]
    print(numpy.mean(lst))

I decided to use PS4 to calculate the running sums of the qscores, then use numpy to find the means. 

I decided I didn't need the init_list function, as that function was just for generating lists correctly in a jupyter notebook.

Had to import gzip and open the data files with it, as they couldn't be unzipped. Also had to add line decoding.
    with gzip.open(file, 'r') as fq:
        line = line.decode('ASCII')

Made a test file with 3 records to see if the histograms were being produced correctly. Had to gzip it to be consistent with the real files.
    gzip test.fq
Now its called test.fq.gz

Ran the script on the test file.

Added argparse arguments to specify read length and file.

PS9 attempt (ended up going back to PS4): 
```
def populate_array(file): 
     """Opens a FASTQ file and decodes Phred quality scores to numbers
     accounting for Phred+33. Sums the quality scores for each position
     and counts the total number of lines in the FASTQ file.
     Returns the array and line count."""
     # qscores = [0] * args.read_len
     with gzip.open(file, 'r') as fq:
         pos = 0
         line_count = 0
         for line in fq:
             line = line.decode('ASCII')
             line = line.strip('\n')
             line_count += 1
             # obtain lines with quality scores
             if line_count % 4 == 0:
                 # specify for position when converting phred score
                 for letter in range(len(line)):
                     qscores[letter, pos] = bioinfo.convert_phred(line[letter])
                 pos +=1
     return (qscores, line_count)
```
I made the script processing.srun to run demux_processing.py on talapas as its a python script, and changed the permissions. 
    chmod 755 processing.srun

I used an interactive node with the command 
    srun --account=bgmp --partition=bgmp --nodes=1 --ntasks-per-node=1 --time=1:00:00 --cpus-per-task=1 --pty bash

The script wasn't running because I didn't add ./ before the file (thought I didn't need to cause it was a python script.)

Once I got the script running, got this error
    qscores = np.zeros((101, line_count//4), dtype=np.int64)
    numpy.core._exceptions._ArrayMemoryError: Unable to allocate 273. GiB for an array with shape (101, 363246735) and data type int64
    Command exited with non-zero status 1

Ended up changing back to the for loop style of PS4 and not using numpy to calucate the means.

processing.srun main command: 
    /usr/bin/time -v ./demux_processing.py  -r 101 -f $file -g $graph

Read 1:
Command being timed: "./demux_processing.py -r 101 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -g 1294_S1_L008_R1_001_dist.png"
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:49:03

Index 1:
Command being timed: "./demux_processing.py -r 8 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -g 1294_S1_L008_R2_001_dist.png"
	Elapsed (wall clock) time (h:mm:ss or m:ss): 20:02.80

Index 2:
Command being timed: "./demux_processing.py -r 8 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -g 1294_S1_L008_R3_001_dist.png"
	Elapsed (wall clock) time (h:mm:ss or m:ss): 21:14.33

Read 2:
Command being timed: "./demux_processing.py -r 101 -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -g 1294_S1_L008_R4_001_dist.png"
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:48:41


Decided on a quality score cut off of 30, as it is typical for selecting high quality data. Also looking at the output histograms, most of the averages were around or above 30, so having a lower threshold likely wouldn't make a major difference in our data. The indexes have only 8 bases, and we discard the read if there's one N, so it makes sense to discard a read that has one low quality value. We are looking at individual qscores rather than the average per read, as higher quality scores may skew the average, and include a lower quality read.

N: amount of indexes that have N's in index files
```
zcat 1294_S1_L008_R2_001.fastq.gz | sed -n 2~4p | grep "N" | wc -l
```
3976613
```
zcat 1294_S1_L008_R3_001.fastq.gz | sed -n 2~4p | grep "N" | wc -l
```
3328051

Total indexes with N's: 3976613 + 3328051 = 7304664

**27 July 2022**

Worked on part 2.

Diagrammed pseudocode options on whiteboard with Lisa, Justine, Jessica, and Kaitlyn.

Wrote out what we needed to record from the algorithm, and the work flow of files.

Decided against a dictionary of known indexes, with the keys as the known indexes and values as the reverse complements of the known indexe

Open all 4 input FASTQ files, read each line by line. 
    for each record, look at lines 1-4
        check if quality score is less than cutoff of 30
            if less than cutoff, write to unknown file
        iterate over the index dict, see if index1 matches a key of the index dict
            if index1 not a key
                write to unknown file
            if index1 is a key
                see if index2 matches a value of the index dict
                    if index2 is not a value
                        write to unknown file
                    if index2 is not a value but not linked to same key index
                        write to unmatched file
                    index2 is value of index1
                        write to matched file



Planned how we might code modifying the header lines to include the indexes.
    for header line
    fh.write(f"{read1[0]} {index1[1]} {index2}\n")
    fh.write(f"{read1[1]}\n{read1[2]}\n{read1[3]}\n") # lines 2-4
    or str: index1 or index2[1]

Planned how it could work saving the lines of the record as an array, I decided against this.
    only change header when writing out to file
    array[0] = array[0] + index1 + rev comp index2             
        
    the arrays will empty each time
    read 1
    [header seq  + Q]
    read 2
    array is only to write to file
    change header line of array to @header_index1_revcompindex2

Decided to use a set rather than a dictionary, as looking up by values is computationally slow. Debated doing two sets, one with the barcodes and one with the reverse complements, or having the forward and reverse complements all in one set, in attempt to reduce the number of times the reverse complement function is called. Realized the reverse complement function must be called anyway to modify the headers with the indexes, so decided on a set with just the forward barcodes, and the reverse complement function will be called for index2

Decided we didn't need to filter for N, as when comparing the barcode to the set of known barcodes, ones with N's would be filtered out and written to the unknown file.

Decided to filter for low quality reads last, as it will be computationally slow and this way it's not trying to convert the quality score for every read.


**28 July 2022**


Note to self:
can use qual score function from bioinfo.py if I decide I want to do average of qscores rather than individual for low qual

Pushed to github. Had to push graphs so the markdown syntax to display them would work.

Note to self:
zip test files
add indexes to headers for test files

**29 July 2022**

Wrote unit tests for high-level functions.

Removed commented out code. 

Wrote input FASTA test files and expected output.

**30 July 2022**

Moved psuedocode from answers.md to psuedocode.txt for easy peer reviewing.

Pushed to github.

**1 August 2022**

Wrote reverse complement function.

    # rev_comp1 = rev_str.replace('A', 'T')
    # print('rv1', rev_comp1)
    # rev_comp2 = rev_comp1.replace('G', 'C')
    # print('rv2', rev_comp2)
    # rev_comp3 = rev_comp2.replace('C', 'G')
    # print('rv3', rev_comp3)
    # rev_comp4 = rev_comp3.replace('T', 'A')
    # print('rv4', rev_comp4)
    
    # for k in rev_str:
    #     if k == 'A':
    #         one = rev_str.replace('A', 'T', 1)
    #     elif k == 'T':
    #         two = rev_str.replace('T', 'A', 1)
    #     elif k == 'C':
    #         three = rev_str.replace('C', 'G', 1)
    #     elif k == 'G':
    #         four = rev_str.replace('G', 'C', 1)
            
    # for k in rev_str:
    #     if k == 'A':
    #         rev_str.replace('A', 'T', 1)
    #     elif k == 'T':
    #         rev_str.replace('T', 'A', 1)
    #     elif k == 'C':
    #         rev_str.replace('C', 'G', 1)
    #     elif k == 'G':
    #         rev_str.replace('G', 'C', 1)

**2 August 2022**

Started working on part 3.

File "/gpfs/projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/demux_alg.py", line 16, in <module>
    valid_indexes = set('GTAGCGTA', 'CGATCGAT', 'GATCAAGG', 'AACAGCGA', 'TAGCCATG', 'CGGTAATC', 'CTCTGGAT', 'TACCGGAT', 'CTAGCTCA',
TypeError: set expected at most 1 argument, got 24

changed to {}

[rferina@talapas-ln1 Assignment-the-third]$ python demux_alg.py -r1 ../TEST-input_FASTQ/r1_test.fq.gz -i1 ../TEST-input_FASTQ/index1_test.fq.gz -i2 ../TEST-input_FASTQ/index2_test.fq.gz  -r2 ../TEST-input_FASTQ/r2_test.fq.gz 
Traceback (most recent call last):
  File "/gpfs/projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/demux_alg.py", line 29, in <module>
    header_r1 = header_r1.strip('\n')
TypeError: a bytes-like object is required, not 'str'

command to run:
 python demux_alg.py -r1 ../TEST-input_FASTQ/r1_test.fq.gz -i1 ../TEST-input_FASTQ/index1_test.fq.gz -i2 ../TEST-input_FASTQ/index2_test.fq.gz  -r2 ../TEST-input_FASTQ/r2_test.fq.gz 

 was initially reading in files and making records lines like this
 # obtain r1 sequence line, strip newline
                        seq_r1 = read1.readline()
                        seq_r1 = seq_r1.strip('\n')
                        # obtain r1 spacer line, strip newline
                        spacer_r1 = read1.readline()
                        spacer_r1 = spacer_r1.strip('\n')
                        # obtain r1 quality score line, strip line
                        quality_r1 = read1.readline()
                        quality_r1 = quality_r1.strip('\n')
                        # print('read1')
                        # print(header_r1)
                        # print(seq_r1)
                        # print(spacer_r1)
                        # print(quality_r1)
                        # obtain i1 sequence line, strip newline
                        seq_i1 = index1.readline()
                        seq_i1 = seq_i1.strip('\n')
                        # obtain i1 spacer line, strip newline
                        spacer_i1 = index1.readline()
                        spacer_i1 = spacer_i1.strip('\n')
                        # obtain i1 quality score line, strip line
                        quality_i1 = index1.readline()
                        quality_i1 = quality_i1.strip('\n')
                        # print('index1')
                        # print(header_i1)
                        # print(seq_i1)
                        # print(spacer_i1)
                        # print(quality_i1)
                        # obtain i2 sequence line, strip newline
                        seq_i2 = index2.readline()
                        seq_i2 = seq_i2.strip('\n')
                        # obtain i2 spacer line, strip newline
                        spacer_i2 = index2.readline()
                        spacer_i2 = spacer_i2.strip('\n')
                        # obtain i2 quality score line, strip line
                        quality_i2 = index2.readline()
                        quality_i2 = quality_i2.strip('\n')
                        # print('index2')
                        # print(header_i2)
                        # print(seq_i2)
                        # print(spacer_i2)
                        # print(quality_i2)
                        # obtain r2 sequence line, strip newline
                        seq_r2 = read2.readline()
                        seq_r2 = seq_r2.strip('\n')
                        # obtain r2 spacer line, strip newline
                        spacer_r2 = read2.readline()
                        spacer_r2 = spacer_r2.strip('\n')
                        # obtain r2 quality score line, strip line
                        quality_r2 = read2.readline()
                        quality_r2 = quality_r2.strip('\n')
                        # print('read2')
                        # print(header_r2)
                        # print(seq_r2)
                        # print(spacer_r2)
                        # print(quality_r2)

learned we could do it in one line with array 

had to add N as a key and value to comp_dict in reverse_complement to avoid keyerror

wrote modify header function

changed TEST-output_FASTQ files to match with the specified header format (- instead of _).

wrote write to output file format, had to remove '\n' after header line bc it already had one from modify header

if doing each score, what if one score 29, want to keep 29; do average of qual score line

tried to write out to unknown test output file, but was using the real valid indexes so everything was getting written to unknown.

realized my test files had incorrect headers, index 2 should be rev comp so changed them.

It worked:
if index1_record[1] in valid_indexes:
                            print('yay')
                        else:
                            read1_indexes_header = modify_header(read1_record[0], index1_record[1], index2_record[1])
                            print(read1_indexes_header)
                            write_out(r1_unknown, read1_indexes_header, read1_record[1], read1_record[2], read1_record[3])

output:
yay
yay
@K00337:83:HJKJNBBXX:8:1101:1347:1191 1:N:0:1 NTCC-NACC

diff -s ../TEST-output_FASTQ/r1_unknown.fq read1_unknown.fq 
Files ../TEST-output_FASTQ/r1_unknown.fq and read1_unknown.fq are identical

however, when checking if both index1 and index2 are not in the valid indexes, my unknown file contains all the sequences.

also my unmatched files are blank, despite my test files containing a high quality unmatched case.

learned should use a dictionary to add the indexes to the output file names.

I was hardcoding the indexes in a set, but need to change the code so that an input index file can be used, and the code can be generalizable.

For these samples, the indexes are in indexex.txt, located in /projects/bgmp/shared/2017_sequencing/.

 python demux_alg.py -r1 ../TEST-input_FASTQ/r1_test.fq.gz -i1 ../TEST-input_FASTQ/index1_test.fq.gz -i2 ../TEST-input_FASTQ/index2_test.fq.gz  -r2 ../TEST-input_FASTQ/r2_test.fq.gz -index /projects/bgmp/shared/2017_sequencing/indexes.txt

 version that doesn't work for unmatched
 else:
                            rev_comp_index2 = bioinfo.reverse_complement(index2_record[1])
                            if rev_comp_index2 in valid_indexes:
                                # check if index1 meets quality score cutoff
                                if bioinfo.qual_score(index1_record[1]) < 30:
                                    print(bioinfo.qual_score(index1_record[1]))
                                    r1_indexes_header = modify_header(read1_record[0], index1_record[1], index2_record[1])
                                    write_out(r1_unknown, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                                    r2_indexes_header = modify_header(read2_record[0], index1_record[1], index2_record[1])
                                    write_out(r2_unknown, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
                                # check if index2 meets quality score cutoff
                                if bioinfo.qual_score(index2_record[1]) < 30:
                                    print(bioinfo.qual_score(index2_record[1]))
                                    r1_indexes_header = modify_header(read1_record[0], index1_record[1], index2_record[1])
                                    write_out(r1_unknown, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                                    r2_indexes_header = modify_header(read2_record[0], index1_record[1], index2_record[1])
                                    write_out(r2_unknown, indexes_header, read2_record[1], read2_record[2], read2_record[3])
                                else:
                                    # print(bioinfo.qual_score(index1_record(1)))
                                    if index1_record[1] == rev_comp_index2:
                                        r1_indexes_header = modify_header(read1_record[0], index1_record[1], index2_record[1])
                                        write_out(r1_match, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                                        r2_indexes_header = modify_header(read2_record[0], index1_record[1], index2_record[1])
                                        write_out(r2_match, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
                                    # if valid indexes don't match 
                                    # elif index1_record[1] != rev_comp_index2:
                                    else:
                                        r1_indexes_header = modify_header(read1_record[0], index1_record[1], index2_record[1])
                                        write_out(r1_unmatch, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                                        r2_indexes_header = modify_header(read2_record[0], index1_record[1], index2_record[1])
                                        write_out(r2_unmatch, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
                        
                                    # checking mathcing after qual score
                                    # if index1_record[1] == rev_comp_index2:
                                    #     r1_indexes_header = modify_header(read1_record[0], index1_record[1], index2_record[1])
                                    #     write_out(r1_match, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                                    #     r2_indexes_header = modify_header(read2_record[0], index1_record[1], index2_record[1])
                                    #     write_out(r2_match, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
                                    # # if valid indexes don't match 
                                    # # elif index1_record[1] != rev_comp_index2:
                                    # else:
                                    #     r1_indexes_header = modify_header(read1_record[0], index1_record[1], index2_record[1])
                                    #     write_out(r1_unmatch, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                                    #     r2_indexes_header = modify_header(read2_record[0], index1_record[1], index2_record[1])
                                    #     write_out(r2_unmatch, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])


changed the order, now checking for valid indexes first and checking for index hopping before checking for the quality score cutoff.

Now my unmatched files are populating, but the index2 part of the header is incorrect.
diff -s r1_notmatch.fq ../TEST-output_FASTQ/r1_unmatched.fq 
1c1
< @K00337:83:HJKJNBBXX:8:1101:1286:1191 1:N:0:1 GCGC-CGCG
---
> @K00337:83:HJKJNBBXX:8:1101:1286:1191 1:N:0:1 GCGC-CTCT


 diff -s r2_notmatch.fq ../TEST-output_FASTQ/r2_unmatched.fq
1c1
< @K00337:83:HJKJNBBXX:8:1101:1286:1191 4:N:0:1 GCGC-CGCG
---
> @K00337:83:HJKJNBBXX:8:1101:1286:1191 4:N:0:1 GCGC-CTCT

The same record is also being written to my unknown files twice.
read2_unknown.fq:
@K00337:83:HJKJNBBXX:8:1101:1347:1191 4:N:0:1 NTCC-NACC
NAAATG
+
#A!JJJ@K00337:83:HJKJNBBXX:8:1101:1347:1191 1:N:0:1 NTCC-NACC
NAAATG
+
#A!JJJ

changed last statement in else to be 
    if index1_record[1] not in valid_indexes or index2_record[1] not in valid_indexes:
rather than two if statements, and now the unknown files don't repeat the same record.

The sequence of r2_unknown didn't match, but I realized my test file was accidentally the same as the r1, so I fixed it. Now they're the same.
 diff -s read2_unknown.fq ../TEST-output_FASTQ/r2_unknown.fq
Files read2_unknown.fq and ../TEST-output_FASTQ/r2_unknown.fq are identical

diff -s read1_unknown.fq ../TEST-output_FASTQ/r1_unknown.fq
Files read1_unknown.fq and ../TEST-output_FASTQ/r1_unknown.fq are identical

diff -s r1_AAAA_match.fq ../TEST-output_FASTQ/AAAA_r1_matched.fq
Files r1_AAAA_match.fq and ../TEST-output_FASTQ/AAAA_r1_matched.fq are identical

diff -s r2_AAAA_match.fq ../TEST-output_FASTQ/AAAA_r2_matched.fq
Files r2_AAAA_match.fq and ../TEST-output_FASTQ/AAAA_r2_matched.fq are identical

troubleshooting by printing the headers, index2, and rev comp of index2, but it isn't actually the rev comp
@K00337:83:HJKJNBBXX:8:1101:1286:1191 1:N:0:1 GCGC-CGCG
CGCG
revcomp CGCG
@K00337:83:HJKJNBBXX:8:1101:1286:1191 4:N:0:1 GCGC-CGCG

CGCG
revcomp2 CGCG
43.0

Something is wrong, index 1 is supposed to be GCGC, but index 2 is supposed to be AGAG, and its reverse comp is CTCT.

Checked reverse_complement of AGAG in bioinfo.py, it works
python bioinfo.py 
CTCT


when printing the index2_record array, the second record is wrong.
['@K00337:83:HJKJNBBXX:8:1101:1265:1191 3:N:0:1' 'TTTT' '+' 'JJJJ']
['@K00337:83:HJKJNBBXX:8:1101:1286:1191 3:N:0:1' 'CGCG' '+' 'JJJJ']
['@K00337:83:HJKJNBBXX:8:1101:1347:1191 3:N:0:1' 'GGTN' '+' '#!!J']

Verified the other record arrays matched the other input files.
read 2:
['@K00337:83:HJKJNBBXX:8:1101:1265:1191 3:N:0:1' 'TTTT' '+' 'JJJJ']
['@K00337:83:HJKJNBBXX:8:1101:1286:1191 3:N:0:1' 'CGCG' '+' 'JJJJ']
['@K00337:83:HJKJNBBXX:8:1101:1347:1191 3:N:0:1' 'GGTN' '+' '#!!J']

my zipped test file was wrong, I didn't have AGAG, so I changed it. Now it passes!

diff -s r1_notmatch.fq ../TEST-output_FASTQ/r1_unmatched.fq
Files r1_notmatch.fq and ../TEST-output_FASTQ/r1_unmatched.fq are identical

diff -s r2_notmatch.fq ../TEST-output_FASTQ/r2_unmatched.fq
Files r2_notmatch.fq and ../TEST-output_FASTQ/r2_unmatched.fq are identical

shouldn't need this code anymore
   # if index2_record[1] not in valid_indexes:
                            #     r1_indexes_header = modify_header(read1_record[0], index1_record[1], index2_record[1])
                            #     write_out(r1_unknown, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                            #     r2_indexes_header = modify_header(read2_record[0], index1_record[1], index2_record[1])
                            #     write_out(r2_unknown, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])   
                        
                                    # checking mathcing after qual score
                                    # if index1_record[1] == rev_comp_index2:
                                    #     r1_indexes_header = modify_header(read1_record[0], index1_record[1], index2_record[1])
                                    #     write_out(r1_match, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                                    #     r2_indexes_header = modify_header(read2_record[0], index1_record[1], index2_record[1])
                                    #     write_out(r2_match, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
                                    # # if valid indexes don't match 
                                    # # elif index1_record[1] != rev_comp_index2:
                                    # else:
                                    #     r1_indexes_header = modify_header(read1_record[0], index1_record[1], index2_record[1])
                                    #     write_out(r1_unmatch, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
                                    #     r2_indexes_header = modify_header(read2_record[0], index1_record[1], index2_record[1])
                                    #     write_out(r2_unmatch, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])


shouldn't need 
# dict comprehension
# {valid_indexes:[open(args.index_r1), open(r2)], for index in indexes_dict}
# files_dict = {valid_indexes:[open(args.index_r1), open(r2)]}

previously had 
write_out(r1_match, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
write_out(r2_match, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
write_out(r1_unknown, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
write_out(r2_unknown, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
write_out(r1_unmatch, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
write_out(r2_unmatch, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])


now am writing out to the files with varying index headings, indexing for read1 or read 2:
write_out(files_dict[index1_record[1]][0], r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])

confirmed this method worked:
$ diff -s AAAA_r1_matched.fq  r1_AAAA_match.fq 
Files AAAA_r1_matched.fq and r1_AAAA_match.fq are identical
diff -s AAAA_r2_matched.fq  r2_AAAA_match.fq 
Files AAAA_r2_matched.fq and r2_AAAA_match.fq are identical

originally had this line when seeing if the key matches the index to open the file, don't need it
# if index1_record[1] in files_dict:

diff -s R1_unknown.fq read1_unknown.fq
Files R1_unknown.fq and read1_unknown.fq are identical

don;t need anymore

r1_unknown = open('read1_unknown.fq', 'w')
r2_unknown = open('read2_unknown.fq', 'w')
r1_match = open('r1_AAAA_match.fq', 'w')
r2_match = open('r2_AAAA_match.fq', 'w')
r1_unmatch = open('r1_notmatch.fq', 'w')
r2_unmatch = open('r2_notmatch.fq', 'w')

# r1_unknown.close()      
# r2_unknown.close()
# r1_match.close()
# r2_match.close()
# r1_unmatch.close()
# r2_unmatch.close()

# valid_indexes = {'GTAGCGTA', 'CGATCGAT', 'GATCAAGG', 'AACAGCGA', 'TAGCCATG', 'CGGTAATC', 'CTCTGGAT', 'TACCGGAT', 'CTAGCTCA',
# 'CACTTCAC', 'GCTACTCT', 'ACGATCAG', 'TATGGCAC', 'TGTTCCGT', 'GTCCTAAG', 'TCGACAAG', 'TCTTCGAC', 'ATCATGCG',
# 'ATCGTGGT',	'TCGAGAGT', 'TCGGATTC', 'GATCTTGC', 'AGAGTCCA', 'AGGATAGC'}


made graph, calculated percentages with a new dict by summing the values of index_count_dict

# how many times this index mismatch with this index
# which unknown, put in bucket unknown
# put results in markdown

/usr/bin/time: cannot run ./demux_alg.py: Permission denied
chmod 755 demux_alg.py
chmod 755 run_demux.sh

  File "/gpfs/projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/./demux_alg.py", line 57, in <module>
    files_dict[index] = [open(f'{args.outputdirectory}/{index}_r1_matched.fq', 'w'), open(f'{args.outputdirectory}/{index}_r2_matched.fq', 'w')]
FileNotFoundError: [Errno 2] No such file or directory: '/projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results%j/GATCAAGG_r1_matched.fq'

removed %j

ran on talapas, 5 hours in noticed my output files were incorrect, as the header isn't written to a newline:
head TACCGGAT_r1_matched.fq
@K00337:83:HJKJNBBXX:8:1101:1174:1701 1:N:0:1 TACCGGAT-TACCGGAT
GCCGACTTNAAATCTGGAAACTTGATGCCTGCAGTCTCAGGGACAGCCGGCTTCTCAGCCTCCTTCTTTTGGGGTAGCTTCTTCTTGGGAGCCTTGCGTTC
+
AAAFFJJJ#JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAJJJJJJFAJJFJJJJJJ<<FA<FFJJJJFJAJJFFJJJJJJJJJJJFAJJAJF77@K00337:83:HJKJNBBXX:8:1101:1296:1701 1:N:0:1 TACCGGAT-TACCGGAT
CTTCCATAAATGTAGGAAATGAAGTCATCAATTTCTATGGTGAACACACAATTCAGCTTGGGGATGACTTTTAGGTGGTCTTCTCTGGTGATGATCTTAAA
+
AAFFFJJJJJJJJJJJJJJJJJJAFAJFJJJJJ-FJJJJJJJJJJJJJJJJJ<FJJJFFJJJJJJJJJJJJJJJJFFFAFJFJJFJJJJJJJJJFJJJJJJ@K00337:83:HJKJNBBXX:8:1101:1336:1701 1:N:0:1 TACCGGAT-TACCGGAT
CTTCTTCTCAAGGTCTTCTCTTTTAGCAGCAAGTGCTTGTTGTTCTGCATTCTCCTCTCTCCACAAGTAATTTCTTTCACTTTGTAGTTCATCTTTCTTAT

output in demux_21828260.out looks good, will also add number of matched records as a sanity check
Number of unmatched records: 23623721
Number of unknown records: 7867981
Matched index counts: {'TAGCCATG': 10629633, 'TCGAGAGT': 11741547, 'CACTTCAC': 4191388, 'AGAGTCCA': 11316780, 'TCTTCGAC': 42094112, 'TGTTCCGT': 15733007, 'GCTACTCT': 7416557, 'CGGTAATC': 5064906, 'GATCAAGG': 6587100, 'TATGGCAC': 11184304, 'TCGACAAG': 3853350, 'CTCTGGAT': 34976387, 'TACCGGAT': 76363857, 'GTCCTAAG': 8830276, 'ATCGTGGT': 6887592, 'AACAGCGA': 8872034, 'CGATCGAT': 5604966, 'TCGGATTC': 4611350, 'CTAGCTCA': 17332036, 'ATCATGCG': 10087503, 'ACGATCAG': 7942853, 'GTAGCGTA': 8119243, 'AGGATAGC': 8673180, 'GATCTTGC': 3641072}
Matched index percentages: {'TAGCCATG': 3.2040608107368187, 'TCGAGAGT': 3.5392219656242565, 'CACTTCAC': 1.263398466663202, 'AGAGTCCA': 3.4111856262328355, 'TCTTCGAC': 12.688311498803998, 'TGTTCCGT': 4.742356689431144, 'GCTACTCT': 2.235552218434618, 'CGGTAATC': 1.5267005760844028, 'GATCAAGG': 1.9855312941100127, 'TATGGCAC': 3.371253752765222, 'TCGACAAG': 1.1615046093362509, 'CTCTGGAT': 10.542835381792083, 'TACCGGAT': 23.018145741288574, 'GTCCTAAG': 2.661685617893851, 'ATCGTGGT': 2.076107764731334, 'AACAGCGA': 2.6742726160841697, 'CGATCGAT': 1.6894893648832752, 'TCGGATTC': 1.3899864482236808, 'CTAGCTCA': 5.224347568526564, 'ATCATGCG': 3.0406480675758125, 'ACGATCAG': 2.394192162866147, 'GTAGCGTA': 2.4473609116278277, 'AGGATAGC': 2.614332606070817, 'GATCTTGC': 1.0975182402131034}
Unmatched: {'TAGCCATG': 389328, 'TCGAGAGT': 473347, 'CACTTCAC': 252184, 'AGAGTCCA': 412398, 'TCTTCGAC': 2291710, 'TGTTCCGT': 570881, 'GCTACTCT': 369326, 'CGGTAATC': 192089, 'GATCAAGG': 268938, 'TATGGCAC': 722204, 'TCGACAAG': 196090, 'CTCTGGAT': 1378876, 'TACCGGAT': 5526459, 'GTCCTAAG': 414410, 'ATCGTGGT': 265830, 'AACAGCGA': 427297, 'CGATCGAT': 187941, 'TCGGATTC': 165990, 'CTAGCTCA': 547172, 'ATCATGCG': 384561, 'ACGATCAG': 255336, 'GTAGCGTA': 361803, 'AGGATAGC': 339343, 'GATCTTGC': 132965}

changed 
r1_indexes_header = modify_header(read1_record[0], index1_record[1], index2_record[1])
to

had to change 
columns = line.split('\t')
to just line.split to work on my test indexes file.

old modify header function calculating the rev comp inside it (less efficient)

def modify_header(line_1, index_1, index_2):
    '''
    Takes in the first line of the record (the header line)
    and index 1 and index 2. Converts index 2 to its reverse
    complement, and adds index 1 and the reverse compliment of 
    index 2 to the end of the line_1 header. Returns the new header.
    '''
    line_1 = str(line_1)
    rev_index_2 = bioinfo.reverse_complement(index_2)
    labeled_header = line_1 + ' ' + index_1 + '-' + rev_index_2 + '\n'
    return labeled_header

read going to unknown rather than unmatched:
python demux_alg.py -r1 ../TEST-input_FASTQ/r1_test.fq.gz -i1 ../TEST-input_FASTQ/index1_test.fq.gz -i2 ../TEST-input_FASTQ/index2_test.fq.gz  -r2 ../TEST-input_FASTQ/r2_test.fq.gz -index ../TEST-input_FASTQ/indexes_for_test.txt -outdir generated_TEST-output/
{'AAAA', 'GCGC', 'ATTC', 'CCCC'}
Number of unmatched records: 1
Number of unknown records: 3
Number of total matched records: 2
Matched index counts: {'AAAA': 1, 'GCGC': 0, 'ATTC': 0, 'CCCC': 1}
Matched index percentages: {'AAAA': 50.0, 'GCGC': 0.0, 'ATTC': 0.0, 'CCCC': 50.0}
Unmatched: {'AAAA': 0, 'GCGC': 1, 'ATTC': 0, 'CCCC': 0}

diff -s ../TEST-output_FASTQ/r1_unknown.fq generated_TEST-output/R1_unknown.fq
8c8
< !#J##J
\ No newline at end of file
---
> !#J##J
(base) [rferina@talapas-ln1 Assignment-the-third]$ diff -s ../TEST-output_FASTQ/r1_unknown.fq generated_TEST-output/R1_unknown.fq
Files ../TEST-output_FASTQ/r1_unknown.fq and generated_TEST-output/R1_unknown.fq are identical

Command being timed: "./demux_alg.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -index /projects/bgmp/shared/2017_sequencing/indexes.txt -outdir /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results"
	User time (seconds): 10709.07
	System time (seconds): 33.55
	Percent of CPU this job got: 94%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:08:30

     JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          21829481      bgmp demux_%j  rferina  R    7:11:51      1 n278

Number of unmatched records: 23623721
Number of unknown records: 7867981
Number of total matched records: 304980270
Matched index counts: {'TCGGATTC': 4163314, 'TCGAGAGT': 10658212, 'AGAGTCCA': 10378366, 'GATCTTGC': 3425453, 'GTCCTAAG': 8164223, 'CACTTCAC': 3833640, 'TCTTCGAC': 39149148, 'TATGGCAC': 10195805, 'AACAGCGA': 8178191, 'CTAGCTCA': 16162895, 'GTAGCGTA': 7450201, 'TCGACAAG': 3548541, 'ACGATCAG': 7441721, 'CGGTAATC': 4498136, 'CTCTGGAT': 32163349, 'TGTTCCGT': 14786868, 'TACCGGAT': 69307073, 'CGATCGAT': 5225776, 'GATCAAGG': 6085915, 'GCTACTCT': 6610857, 'AGGATAGC': 8078057, 'ATCATGCG': 9264615, 'ATCGTGGT': 6357656, 'TAGCCATG': 9852258}
Matched index percentages: {'TCGGATTC': 1.3651092905124649, 'TCGAGAGT': 3.4947218061024077, 'AGAGTCCA': 3.402963083480777, 'GATCTTGC': 1.123172000601875, 'GTCCTAAG': 2.6769675953136245, 'CACTTCAC': 1.2570124618225302, 'TCTTCGAC': 12.83661661129751, 'TATGGCAC': 3.343103145655947, 'AACAGCGA': 2.6815475637161708, 'CTAGCTCA': 5.299652662777169, 'GTAGCGTA': 2.4428468766192646, 'TCGACAAG': 1.1635313326989973, 'ACGATCAG': 2.440066368883469, 'CGGTAATC': 1.4748940972476678, 'CTCTGGAT': 10.546042535800758, 'TGTTCCGT': 4.848467082805062, 'TACCGGAT': 22.725100545028702, 'CGATCGAT': 1.71348002282246, 'GATCAAGG': 1.9955110538789935, 'GCTACTCT': 2.1676343194266305, 'AGGATAGC': 2.648714620129361, 'ATCATGCG': 3.0377751977201672, 'ATCGTGGT': 2.084612227538522, 'TAGCCATG': 3.2304574981194683}
Unmatched: {'TCGGATTC': 282016, 'TCGAGAGT': 870049, 'AGAGTCCA': 715130, 'GATCTTGC': 219113, 'GTCCTAAG': 580802, 'CACTTCAC': 336290, 'TCTTCGAC': 3260293, 'TATGGCAC': 1074097, 'AACAGCGA': 585239, 'CTAGCTCA': 913700, 'GTAGCGTA': 508995, 'TCGACAAG': 290019, 'ACGATCAG': 492084, 'CGGTAATC': 294010, 'CTCTGGAT': 2129147, 'TGTTCCGT': 1067764, 'TACCGGAT': 7248375, 'CGATCGAT': 317203, 'GATCAAGG': 400298, 'GCTACTCT': 532220, 'AGGATAGC': 526332, 'ATCATGCG': 631508, 'ATCGTGGT': 462526, 'TAGCCATG': 594251}



test with matching graph, rotated axis label
python demux_alg.py -r1 ../TEST-input_FASTQ/r1_test.fq.gz -i1 ../TEST-input_FASTQ/index1_test.fq.gz -i2 ../TEST-input_FASTQ/index2_test.fq.gz  -r2 ../TEST-input_FASTQ/r2_test.fq.gz -index ../TEST-input_FASTQ/indexes_for_test.txt -outdir generated_TEST-output/
Number of unmatched records: 2
Number of unknown records: 2
Number of total matched records: 2
Matched index counts: {'GCGC': 0, 'CCCC': 1, 'AAAA': 1, 'GAAT': 0}
Matched index percentages: {'GCGC': 0.0, 'CCCC': 50.0, 'AAAA': 50.0, 'GAAT': 0.0}
Unmatched: {'GCGC': 1, 'CCCC': 0, 'AAAA': 0, 'GAAT': 1}