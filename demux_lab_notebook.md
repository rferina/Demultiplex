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
Can use qual_score function from bioinfo.py if I decide I want to do average of qscores rather than individual for low quality

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

changed to {} instead of set() and it worked.

Ran demux_alg.py on my input test files with this command:
 python demux_alg.py -r1 ../TEST-input_FASTQ/r1_test.fq.gz -i1 ../TEST-input_FASTQ/index1_test.fq.gz -i2 ../TEST-input_FASTQ/index2_test.fq.gz  -r2 ../TEST-input_FASTQ/r2_test.fq.gz 

 was initially reading in files and making record lines like this
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

It is much more efficient to code with an array to save the lines of a record, so changed to an array after hearing my peer's opinions.  

Realized I had to add N as a key and value to my comp_dict in my reverse_complement function to avoid a keyerror.

Wrote a function modify_header to add index 1 and the reverse complement of index 2 to the header line. Realized the format is supposed to be index1-revcompindex2, not index1_revcompindex2, so I changed TEST-output_FASTQ files to match with the specified header format.

Wrote write_out function to write out the read to an output file, had to remove '\n' after the first line because it's the header line and it already had one from modify_header.

Debated using the average of qscores for the cutoff, or individual qscores in the indexes. If using individual scores, what if one score was 29 and the rest were high? We wouldn't want to discard that read as the rest is high quality. So I decided to change to the average of the quality score line needing to meet the cutoff, which gives some leeway for an otherwise high quality read.

Tried to write out to the test output files, but was using the real valid indexes so everything was getting written to unknown. Made a set of my valid test indexes.

Realized my test files had incorrect headers, as index 2 should be the rev complement so I changed them.

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

Checked that my output files matched the output files I defined.
diff -s ../TEST-output_FASTQ/r1_unknown.fq read1_unknown.fq 
Files ../TEST-output_FASTQ/r1_unknown.fq and read1_unknown.fq are identical

However, when checking if both index1 and index2 are not in the valid indexes, my unknown file contains all the sequences.

Also my unmatched files are blank, despite my test files containing a high quality unmatched case.

After discussing with my peers, I learned I should use a dictionary to add the indexes to the output file names.

I was hardcoding the indexes in a set, but need to change the code so that an input index file can be used, and the code can be generalizable.

**3 August 2022**

For these real samples, learned that the indexes are in index.txt, located in /projects/bgmp/shared/2017_sequencing/. Ended up adding an argparse arugment and reading in the file instead of hardcoding the index set.
 python demux_alg.py -r1 ../TEST-input_FASTQ/r1_test.fq.gz -i1 ../TEST-input_FASTQ/index1_test.fq.gz -i2 ../TEST-input_FASTQ/index2_test.fq.gz  -r2 ../TEST-input_FASTQ/r2_test.fq.gz -index /projects/bgmp/shared/2017_sequencing/indexes.txt


 Version that doesn't work for unmatched
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


Changed the order, now checking for valid indexes first and checking for index hopping before checking for the quality score cutoff.

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

Changed last statement in the else to be: 
    if index1_record[1] not in valid_indexes or index2_record[1] not in valid_indexes:
rather than two if statements, and now the unknown files don't repeat the same record.

The sequence of r2_unknown didn't match, but I realized my test file was accidentally the same as the r1, so I fixed it. Now they match.
 diff -s read2_unknown.fq ../TEST-output_FASTQ/r2_unknown.fq
Files read2_unknown.fq and ../TEST-output_FASTQ/r2_unknown.fq are identical

diff -s read1_unknown.fq ../TEST-output_FASTQ/r1_unknown.fq
Files read1_unknown.fq and ../TEST-output_FASTQ/r1_unknown.fq are identical

diff -s r1_AAAA_match.fq ../TEST-output_FASTQ/AAAA_r1_matched.fq
Files r1_AAAA_match.fq and ../TEST-output_FASTQ/AAAA_r1_matched.fq are identical

diff -s r2_AAAA_match.fq ../TEST-output_FASTQ/AAAA_r2_matched.fq
Files r2_AAAA_match.fq and ../TEST-output_FASTQ/AAAA_r2_matched.fq are identical

Troubleshooting by printing the headers, index2, and rev comp of index2, but it isn't actually the reverse complement.
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

Learned my zipped test file was wrong, I didn't have AGAG, so I changed it. Now it passes!

diff -s r1_notmatch.fq ../TEST-output_FASTQ/r1_unmatched.fq
Files r1_notmatch.fq and ../TEST-output_FASTQ/r1_unmatched.fq are identical

diff -s r2_notmatch.fq ../TEST-output_FASTQ/r2_unmatched.fq
Files r2_notmatch.fq and ../TEST-output_FASTQ/r2_unmatched.fq are identical

this code is a draft, not needed anymore
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



Previously had 
write_out(r1_match, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
write_out(r2_match, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
write_out(r1_unknown, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
write_out(r2_unknown, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])
write_out(r1_unmatch, r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])
write_out(r2_unmatch, r2_indexes_header, read2_record[1], read2_record[2], read2_record[3])

Mow am writing out to the files with varying index headings, indexing for read1 or read 2:
write_out(files_dict[index1_record[1]][0], r1_indexes_header, read1_record[1], read1_record[2], read1_record[3])

Confirmed this method worked:
$ diff -s AAAA_r1_matched.fq  r1_AAAA_match.fq 
Files AAAA_r1_matched.fq and r1_AAAA_match.fq are identical
diff -s AAAA_r2_matched.fq  r2_AAAA_match.fq 
Files AAAA_r2_matched.fq and r2_AAAA_match.fq are identical

Originally had this line when seeing if the key matches the index to open the file, don't need it
# if index1_record[1] in files_dict:

diff -s R1_unknown.fq read1_unknown.fq
Files R1_unknown.fq and read1_unknown.fq are identical

Don't need this hardcoding opening/closing the files anymore:
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


Made distribution graphs, and calculated percentages with a new dict by summing the values of index_count_dict.

Tried to run with the script run_demux.sh, got this error:
/usr/bin/time: cannot run ./demux_alg.py: Permission denied
Fixed with:
chmod 755 demux_alg.py
chmod 755 run_demux.sh

Tried to make the results file be different for each run but it didn't work, had to remove %j
  File "/gpfs/projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/./demux_alg.py", line 57, in <module>
    files_dict[index] = [open(f'{args.outputdirectory}/{index}_r1_matched.fq', 'w'), open(f'{args.outputdirectory}/{index}_r2_matched.fq', 'w')]
FileNotFoundError: [Errno 2] No such file or directory: '/projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results%j/GATCAAGG_r1_matched.fq'


Ran on talapas, 5 hours in I noticed my output file formatting was incorrect, as the header isn't written to a newline:
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
r1_indexes_header = modify_header(read1_record[0], index1_record[1], rev_comp_index2)
because I already calculated the rev_comp_index2, so I also changed modify_header to input the rev_comp_index2 rather than also converting index2 inside the function unnecessarily.

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


Had to change 
columns = line.split('\t')
to just line.split() to work on my test indexes file.

**4 August 2022**

In test files, one read was going to unknown rather than unmatched:
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

Had to add newline.
(base) [rferina@talapas-ln1 Assignment-the-third]$ diff -s ../TEST-output_FASTQ/r1_unknown.fq generated_TEST-output/R1_unknown.fq
Files ../TEST-output_FASTQ/r1_unknown.fq and generated_TEST-output/R1_unknown.fq are identical

Now passing test files, went back to running on read files:
Command being timed: "./demux_alg.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -index /projects/bgmp/shared/2017_sequencing/indexes.txt -outdir /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results"
	User time (seconds): 10709.07
	System time (seconds): 33.55
	Percent of CPU this job got: 94%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:08:30

     JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
          21829481      bgmp demux_%j  rferina  R    7:11:51      1 n278

Output of run:
Number of unmatched records: 23623721
Number of unknown records: 7867981
Number of total matched records: 304980270
Matched index counts: {'TCGGATTC': 4163314, 'TCGAGAGT': 10658212, 'AGAGTCCA': 10378366, 'GATCTTGC': 3425453, 'GTCCTAAG': 8164223, 'CACTTCAC': 3833640, 'TCTTCGAC': 39149148, 'TATGGCAC': 10195805, 'AACAGCGA': 8178191, 'CTAGCTCA': 16162895, 'GTAGCGTA': 7450201, 'TCGACAAG': 3548541, 'ACGATCAG': 7441721, 'CGGTAATC': 4498136, 'CTCTGGAT': 32163349, 'TGTTCCGT': 14786868, 'TACCGGAT': 69307073, 'CGATCGAT': 5225776, 'GATCAAGG': 6085915, 'GCTACTCT': 6610857, 'AGGATAGC': 8078057, 'ATCATGCG': 9264615, 'ATCGTGGT': 6357656, 'TAGCCATG': 9852258}
Matched index percentages: {'TCGGATTC': 1.3651092905124649, 'TCGAGAGT': 3.4947218061024077, 'AGAGTCCA': 3.402963083480777, 'GATCTTGC': 1.123172000601875, 'GTCCTAAG': 2.6769675953136245, 'CACTTCAC': 1.2570124618225302, 'TCTTCGAC': 12.83661661129751, 'TATGGCAC': 3.343103145655947, 'AACAGCGA': 2.6815475637161708, 'CTAGCTCA': 5.299652662777169, 'GTAGCGTA': 2.4428468766192646, 'TCGACAAG': 1.1635313326989973, 'ACGATCAG': 2.440066368883469, 'CGGTAATC': 1.4748940972476678, 'CTCTGGAT': 10.546042535800758, 'TGTTCCGT': 4.848467082805062, 'TACCGGAT': 22.725100545028702, 'CGATCGAT': 1.71348002282246, 'GATCAAGG': 1.9955110538789935, 'GCTACTCT': 2.1676343194266305, 'AGGATAGC': 2.648714620129361, 'ATCATGCG': 3.0377751977201672, 'ATCGTGGT': 2.084612227538522, 'TAGCCATG': 3.2304574981194683}
Unmatched: {'TCGGATTC': 282016, 'TCGAGAGT': 870049, 'AGAGTCCA': 715130, 'GATCTTGC': 219113, 'GTCCTAAG': 580802, 'CACTTCAC': 336290, 'TCTTCGAC': 3260293, 'TATGGCAC': 1074097, 'AACAGCGA': 585239, 'CTAGCTCA': 913700, 'GTAGCGTA': 508995, 'TCGACAAG': 290019, 'ACGATCAG': 492084, 'CGGTAATC': 294010, 'CTCTGGAT': 2129147, 'TGTTCCGT': 1067764, 'TACCGGAT': 7248375, 'CGATCGAT': 317203, 'GATCAAGG': 400298, 'GCTACTCT': 532220, 'AGGATAGC': 526332, 'ATCATGCG': 631508, 'ATCGTGGT': 462526, 'TAGCCATG': 594251}


Realized I should add more graphs to interpret the statistics, so I tested making a frequency graph and rotating the axis labels to be able to see the indexes on the x axis, using the test input files.
python demux_alg.py -r1 ../TEST-input_FASTQ/r1_test.fq.gz -i1 ../TEST-input_FASTQ/index1_test.fq.gz -i2 ../TEST-input_FASTQ/index2_test.fq.gz  -r2 ../TEST-input_FASTQ/r2_test.fq.gz -index ../TEST-input_FASTQ/indexes_for_test.txt -outdir generated_TEST-output/
Number of unmatched records: 2
Number of unknown records: 2
Number of total matched records: 2
Matched index counts: {'GCGC': 0, 'CCCC': 1, 'AAAA': 1, 'GAAT': 0}
Matched index percentages: {'GCGC': 0.0, 'CCCC': 50.0, 'AAAA': 50.0, 'GAAT': 0.0}
Unmatched: {'GCGC': 1, 'CCCC': 0, 'AAAA': 0, 'GAAT': 1}

Added a command in run_demux.sh to zip the output files, it ran successfully but took 4 hours to zip using gzip.
Command being timed: "./demux_alg.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -index /projects/bgmp/shared/2017_sequencing/indexes.txt -outdir /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results"
	User time (seconds): 11038.44
	System time (seconds): 45.24
	Percent of CPU this job got: 95%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:13:12
	Average shared text size (kbytes): 0

Command being timed: "gzip /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AACAGCGA_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AACAGCGA_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ACGATCAG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ACGATCAG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AGAGTCCA_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AGAGTCCA_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AGGATAGC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AGGATAGC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ATCATGCG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ATCATGCG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ATCGTGGT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ATCGTGGT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CACTTCAC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CACTTCAC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CGATCGAT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CGATCGAT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CGGTAATC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CGGTAATC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CTAGCTCA_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CTAGCTCA_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CTCTGGAT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CTCTGGAT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GATCAAGG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GATCAAGG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GATCTTGC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GATCTTGC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GCTACTCT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GCTACTCT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GTAGCGTA_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GTAGCGTA_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GTCCTAAG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GTCCTAAG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/r1_unknown.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/r1_unmatched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/r2_unknown.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/r2_unmatched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TACCGGAT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TACCGGAT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TAGCCATG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TAGCCATG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TATGGCAC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TATGGCAC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGACAAG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGACAAG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGAGAGT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGAGAGT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGGATTC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGGATTC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCTTCGAC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCTTCGAC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TGTTCCGT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TGTTCCGT_r2_matched.fq"
	User time (seconds): 17082.62
	System time (seconds): 49.48
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 4:46:09


Output of run:
Number of unmatched records: 23623721
Number of unknown records: 7867981
Number of total matched records: 304980270
Matched index counts: {'TCGAGAGT': 10658212, 'CGATCGAT': 5225776, 'GCTACTCT': 6610857, 'AGAGTCCA': 10378366, 'TACCGGAT': 69307073, 'ATCATGCG': 9264615, 'GTAGCGTA': 7450201, 'TGTTCCGT': 14786868, 'AGGATAGC': 8078057, 'ATCGTGGT': 6357656, 'GATCAAGG': 6085915, 'TCGACAAG': 3548541, 'TATGGCAC': 10195805, 'CTAGCTCA': 16162895, 'AACAGCGA': 8178191, 'CTCTGGAT': 32163349, 'TAGCCATG': 9852258, 'CACTTCAC': 3833640, 'TCTTCGAC': 39149148, 'TCGGATTC': 4163314, 'ACGATCAG': 7441721, 'CGGTAATC': 4498136, 'GATCTTGC': 3425453, 'GTCCTAAG': 8164223}
Matched index percentages: {'TCGAGAGT': 3.4947218061024077, 'CGATCGAT': 1.71348002282246, 'GCTACTCT': 2.1676343194266305, 'AGAGTCCA': 3.402963083480777, 'TACCGGAT': 22.725100545028702, 'ATCATGCG': 3.0377751977201672, 'GTAGCGTA': 2.4428468766192646, 'TGTTCCGT': 4.848467082805062, 'AGGATAGC': 2.648714620129361, 'ATCGTGGT': 2.084612227538522, 'GATCAAGG': 1.9955110538789935, 'TCGACAAG': 1.1635313326989973, 'TATGGCAC': 3.343103145655947, 'CTAGCTCA': 5.299652662777169, 'AACAGCGA': 2.6815475637161708, 'CTCTGGAT': 10.546042535800758, 'TAGCCATG': 3.2304574981194683, 'CACTTCAC': 1.2570124618225302, 'TCTTCGAC': 12.83661661129751, 'TCGGATTC': 1.3651092905124649, 'ACGATCAG': 2.440066368883469, 'CGGTAATC': 1.4748940972476678, 'GATCTTGC': 1.123172000601875, 'GTCCTAAG': 2.6769675953136245}
Unmatched: {'TCGAGAGT': 870049, 'CGATCGAT': 317203, 'GCTACTCT': 532220, 'AGAGTCCA': 715130, 'TACCGGAT': 7248375, 'ATCATGCG': 631508, 'GTAGCGTA': 508995, 'TGTTCCGT': 1067764, 'AGGATAGC': 526332, 'ATCGTGGT': 462526, 'GATCAAGG': 400298, 'TCGACAAG': 290019, 'TATGGCAC': 1074097, 'CTAGCTCA': 913700, 'AACAGCGA': 585239, 'CTCTGGAT': 2129147, 'TAGCCATG': 594251, 'CACTTCAC': 336290, 'TCTTCGAC': 3260293, 'TCGGATTC': 282016, 'ACGATCAG': 492084, 'CGGTAATC': 294010, 'GATCTTGC': 219113, 'GTCCTAAG': 580802}

Comparing my output to peers, I realized I had the definition of unknown wrong--if one barcode doesn't match, the whole read should be in the unknown file. Before, I was counting records with one valid index as unmatched instead of unknown.

I also learned I can save time gzipping by using pigz.

10 August 2022

Removed gzip command from run_demux.sh, replaced with pigz.
/usr/bin/time -v gzip $results/*.fq
/usr/bin/time -v pigz $results/*.fq

Made more test cases for my updated definition of unknown reads in my test files.

I previously was calculating the percentage of matched indexes over the total matched indexes, however this is not very informative for demultiplexing. I am changing the denominator to be over the total number of reads, so we can see how many matched out of all the reads.
changed from:
total_matched = sum(index_count_dict.values())
index_percentages = {}
for index in index_count_dict: 
    index_percentages[index] = index_count_dict[index] / total_matched * 100
to: 
#find the total number of matched indexes 
total_matched = sum(index_count_dict.values())
#calculate percentages of matched indexes over all reads
for index in index_count_dict: 
    index_percentages[index] = index_count_dict[index] / total_count * 100


Now that the unknown and percentage are corrected, I reran on talapas. However, it still took 4.5 hours to zip, even with pigz.

Number of unmatched records: 517612
Number of unknown records: 57748853
Number of total matched records: 304980270
Unmatched index counts: {'GATCAAGG': 28066, 'CGGTAATC': 16157, 'TATGGCAC': 192812, 'TAGCCATG': 20535, 'TCGGATTC': 12198, 'GTAGCGTA': 17615, 'TCGACAAG': 27800, 'AACAGCGA': 17929, 'GCTACTCT': 17199, 'AGAGTCCA': 19691, 'CGATCGAT': 14019, 'TACCGGAT': 111390, 'AGGATAGC': 17674, 'GTCCTAAG': 22282, 'ATCATGCG': 26882, 'ACGATCAG': 18369, 'CACTTCAC': 11945, 'ATCGTGGT': 17348, 'CTCTGGAT': 66503, 'CTAGCTCA': 46507, 'TGTTCCGT': 189176, 'TCGAGAGT': 24501, 'TCTTCGAC': 89525, 'GATCTTGC': 9101}
Matched index counts: {'GATCAAGG': 6085915, 'CGGTAATC': 4498136, 'TATGGCAC': 10195805, 'TAGCCATG': 9852258, 'TCGGATTC': 4163314, 'GTAGCGTA': 7450201, 'TCGACAAG': 3548541, 'AACAGCGA': 8178191, 'GCTACTCT': 6610857, 'AGAGTCCA': 10378366, 'CGATCGAT': 5225776, 'TACCGGAT': 69307073, 'AGGATAGC': 8078057, 'GTCCTAAG': 8164223, 'ATCATGCG': 9264615, 'ACGATCAG': 7441721, 'CACTTCAC': 3833640, 'ATCGTGGT': 6357656, 'CTCTGGAT': 32163349, 'CTAGCTCA': 16162895, 'TGTTCCGT': 14786868, 'TCGAGAGT': 10658212, 'TCTTCGAC': 39149148, 'GATCTTGC': 3425453}
Matched index percentages: {'GATCAAGG': 1.6754218038601227, 'CGGTAATC': 1.2383142273804608, 'TATGGCAC': 2.8068538592645575, 'TAGCCATG': 2.7122770972738404, 'TCGGATTC': 1.146139414026667, 'GTAGCGTA': 2.0510028810031837, 'TCGACAAG': 0.9768954977668278, 'AACAGCGA': 2.2514148681886983, 'GCTACTCT': 1.8199356974261585, 'AGAGTCCA': 2.85711198477806, 'CGATCGAT': 1.4386298613255257, 'TACCGGAT': 19.079888770369816, 'AGGATAGC': 2.2238484813910304, 'GTCCTAAG': 2.2475695480098397, 'ATCATGCG': 2.550501933623712, 'ACGATCAG': 2.0486683796345755, 'CACTTCAC': 1.055381819192401, 'ATCGTGGT': 1.7502307350401924, 'CTCTGGAT': 8.854408285321545, 'CTAGCTCA': 4.449563737992029, 'TGTTCCGT': 4.070750422574342, 'TCGAGAGT': 2.9341521817119705, 'TCTTCGAC': 10.777563630406755, 'GATCTTGC': 0.9430099901655}


Command being timed: "./demux_alg.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -index /projects/bgmp/shared/2017_sequencing/indexes.txt -outdir /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results"
	User time (seconds): 10912.74
	System time (seconds): 34.15
	Percent of CPU this job got: 95%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:10:23


Command being timed: "pigz /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AACAGCGA_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AACAGCGA_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ACGATCAG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ACGATCAG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AGAGTCCA_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AGAGTCCA_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AGGATAGC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AGGATAGC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ATCATGCG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ATCATGCG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ATCGTGGT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ATCGTGGT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CACTTCAC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CACTTCAC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CGATCGAT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CGATCGAT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CGGTAATC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CGGTAATC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CTAGCTCA_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CTAGCTCA_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CTCTGGAT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CTCTGGAT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GATCAAGG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GATCAAGG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GATCTTGC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GATCTTGC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GCTACTCT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GCTACTCT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GTAGCGTA_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GTAGCGTA_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GTCCTAAG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GTCCTAAG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/r1_unknown.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/r1_unmatched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/r2_unknown.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/r2_unmatched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TACCGGAT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TACCGGAT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TAGCCATG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TAGCCATG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TATGGCAC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TATGGCAC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGACAAG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGACAAG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGAGAGT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGAGAGT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGGATTC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGGATTC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCTTCGAC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCTTCGAC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TGTTCCGT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TGTTCCGT_r2_matched.fq"
	User time (seconds): 17043.40
	System time (seconds): 21.66
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 4:44:27


**11  August 2022**

I learned I needed to change the number of cpus in order to make pgiz run fast, so I added this line to run_demux.sh:
cpus-per-task=12

I also decided to report the total reads in the file.
#find the total reads in the file
total_count = total_matched + unmatched_count + unknown_count
index_percentages = {}

I realized a dictionary of percentages was not the most informative, so I also added a pie chart. I ran into issues with the barcodes overlapping when displayed, so moved the title of the plot up with the parameter pad=28. I also decided to rotate the barcode labels so they wouldn't overlap with each other with rotatelabels=True, and reported the percentage for each slice with autopct='%1.1f%%'.

I reran on talapas, and my output pie chart slightly overlaps with the title, but I already moved it as high as it could be without cutting it off. Pigz only took 30 minutes with the additional cpus!

Command being timed: "./demux_alg.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -index /projects/bgmp/shared/2017_sequencing/indexes.txt -outdir /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results"
	User time (seconds): 13386.73
	System time (seconds): 73.43
	Percent of CPU this job got: 93%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 4:00:27

Command being timed: "pigz /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AACAGCGA_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AACAGCGA_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ACGATCAG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ACGATCAG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AGAGTCCA_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AGAGTCCA_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AGGATAGC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/AGGATAGC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ATCATGCG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ATCATGCG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ATCGTGGT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/ATCGTGGT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CACTTCAC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CACTTCAC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CGATCGAT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CGATCGAT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CGGTAATC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CGGTAATC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CTAGCTCA_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CTAGCTCA_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CTCTGGAT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/CTCTGGAT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GATCAAGG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GATCAAGG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GATCTTGC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GATCTTGC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GCTACTCT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GCTACTCT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GTAGCGTA_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GTAGCGTA_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GTCCTAAG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/GTCCTAAG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/r1_unknown.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/r1_unmatched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/r2_unknown.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/r2_unmatched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TACCGGAT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TACCGGAT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TAGCCATG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TAGCCATG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TATGGCAC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TATGGCAC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGACAAG_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGACAAG_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGAGAGT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGAGAGT_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGGATTC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCGGATTC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCTTCGAC_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TCTTCGAC_r2_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TGTTCCGT_r1_matched.fq /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results/TGTTCCGT_r2_matched.fq"
	User time (seconds): 25740.05
	System time (seconds): 34.25
	Percent of CPU this job got: 1192%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 36:01.21

Output:
Number of unmatched records: 517612
Number of unknown records: 57748853
Number of total matched records: 304980270
Total Number of Records: 363246735
Unmatched index counts: {'AGGATAGC': 17674, 'AACAGCGA': 17929, 'ACGATCAG': 18369, 'ATCGTGGT': 17348, 'TAGCCATG': 20535, 'CTCTGGAT': 66503, 'AGAGTCCA': 19691, 'TCGAGAGT': 24501, 'CGGTAATC': 16157, 'GATCAAGG': 28066, 'GTCCTAAG': 22282, 'TCGACAAG': 27800, 'ATCATGCG': 26882, 'GTAGCGTA': 17615, 'TGTTCCGT': 189176, 'TCTTCGAC': 89525, 'CTAGCTCA': 46507, 'GATCTTGC': 9101, 'CGATCGAT': 14019, 'TACCGGAT': 111390, 'GCTACTCT': 17199, 'TCGGATTC': 12198, 'TATGGCAC': 192812, 'CACTTCAC': 11945}
Matched index counts: {'AGGATAGC': 8078057, 'AACAGCGA': 8178191, 'ACGATCAG': 7441721, 'ATCGTGGT': 6357656, 'TAGCCATG': 9852258, 'CTCTGGAT': 32163349, 'AGAGTCCA': 10378366, 'TCGAGAGT': 10658212, 'CGGTAATC': 4498136, 'GATCAAGG': 6085915, 'GTCCTAAG': 8164223, 'TCGACAAG': 3548541, 'ATCATGCG': 9264615, 'GTAGCGTA': 7450201, 'TGTTCCGT': 14786868, 'TCTTCGAC': 39149148, 'CTAGCTCA': 16162895, 'GATCTTGC': 3425453, 'CGATCGAT': 5225776, 'TACCGGAT': 69307073, 'GCTACTCT': 6610857, 'TCGGATTC': 4163314, 'TATGGCAC': 10195805, 'CACTTCAC': 3833640}
Matched index percentages: {'AGGATAGC': 2.2238484813910304, 'AACAGCGA': 2.2514148681886983, 'ACGATCAG': 2.0486683796345755, 'ATCGTGGT': 1.7502307350401924, 'TAGCCATG': 2.7122770972738404, 'CTCTGGAT': 8.854408285321545, 'AGAGTCCA': 2.85711198477806, 'TCGAGAGT': 2.9341521817119705, 'CGGTAATC': 1.2383142273804608, 'GATCAAGG': 1.6754218038601227, 'GTCCTAAG': 2.2475695480098397, 'TCGACAAG': 0.9768954977668278, 'ATCATGCG': 2.550501933623712, 'GTAGCGTA': 2.0510028810031837, 'TGTTCCGT': 4.070750422574342, 'TCTTCGAC': 10.777563630406755, 'CTAGCTCA': 4.449563737992029, 'GATCTTGC': 0.9430099901655, 'CGATCGAT': 1.4386298613255257, 'TACCGGAT': 19.079888770369816, 'GCTACTCT': 1.8199356974261585, 'TCGGATTC': 1.146139414026667, 'TATGGCAC': 2.8068538592645575, 'CACTTCAC': 1.055381819192401}


Added the output to summary_stats.md, and added the graphs inline as well.

Added type annotations to my functions.

**12 August 2022**

Realized I was double counting the unmatched indexes because of my old definiton.
changed this:
#count times an index has a mismatch for index1
 if index1_record[1] in unmatched_index_dict:
    unmatched_index_dict[index1_record[1]] += 1
#count times an index has a mismatch for index2
if rev_comp_index2 in unmatched_index_dict:
    unmatched_index_dict[rev_comp_index2] += 1
to:
#count how many mismatches
if str(index1_record[1] + rev_comp_index2) in unmatched_index_dict:
    unmatched_index_dict[str(index1_record[1]) + ':' + rev_comp_index2] += 1
else:
    unmatched_index_dict[str(index1_record[1]) + ':' + rev_comp_index2] = 1

Adding the else statements means I didn't need this initialization.
#count the frequencies of matched indexes
index_count_dict = {}
for index in valid_indexes:
    index_count_dict[index] = 0

Reran on talapas: 
Command being timed: "./demux_alg.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -index /projects/bgmp/shared/2017_sequencing/indexes.txt -outdir /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results"
	User time (seconds): 12777.88
	System time (seconds): 50.97
	Percent of CPU this job got: 93%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:48:07


Output:
Number of unmatched records: 517612
Number of unknown records: 57748853
Number of total matched records: 304980270
Total Number of Records: 363246735
Unmatched index counts: {'TATGGCAC:TGTTCCGT': 1, 'CTAGCTCA:CGGTAATC': 1, 'CTCTGGAT:AGGATAGC': 1, 'TCGACAAG:ATCATGCG': 1, 'TGTTCCGT:TATGGCAC': 1, 'TACCGGAT:TCGAGAGT': 1, 'GATCAAGG:TCTTCGAC': 1, 'ATCATGCG:TATGGCAC': 1, 'GTCCTAAG:TATGGCAC': 1, 'TCTTCGAC:TCGGATTC': 1, 'TACCGGAT:AGGATAGC': 1, 'GATCAAGG:GCTACTCT': 1, 'CACTTCAC:TCTTCGAC': 1, 'TACCGGAT:TCGACAAG': 1, 'CTCTGGAT:TACCGGAT': 1, 'AGAGTCCA:TACCGGAT': 1, 'CACTTCAC:ACGATCAG': 1, 'TATGGCAC:AACAGCGA': 1, 'CGATCGAT:GTAGCGTA': 1, 'TCTTCGAC:CTCTGGAT': 1, 'CTAGCTCA:CTCTGGAT': 1, 'CGATCGAT:TGTTCCGT': 1, 'AGGATAGC:TCTTCGAC': 1, 'TGTTCCGT:TACCGGAT': 1, 'CTAGCTCA:TCGACAAG': 1, 'TACCGGAT:CTCTGGAT': 1, 'TATGGCAC:TCTTCGAC': 1, 'TCGGATTC:ACGATCAG': 1, 'TACCGGAT:TCTTCGAC': 1, 'CGGTAATC:TACCGGAT': 1, 'GATCAAGG:ACGATCAG': 1, 'CGATCGAT:TACCGGAT': 1, 'TCTTCGAC:ATCGTGGT': 1, 'GTCCTAAG:CTAGCTCA': 1, 'GCTACTCT:GTAGCGTA': 1, 'GCTACTCT:TACCGGAT': 1, 'TGTTCCGT:TAGCCATG': 1, 'AGAGTCCA:CGATCGAT': 1, 'TACCGGAT:GCTACTCT': 1, 'TCGACAAG:AGAGTCCA': 1, 'AGGATAGC:CTCTGGAT': 1, 'GATCAAGG:GTCCTAAG': 1, 'ATCATGCG:TCTTCGAC': 1, 'AGGATAGC:TACCGGAT': 1, 'GTCCTAAG:TACCGGAT': 1, 'AGGATAGC:GTAGCGTA': 1, 'AGGATAGC:GCTACTCT': 1, 'CTCTGGAT:AACAGCGA': 1, 'CTCTGGAT:GCTACTCT': 1, 'CACTTCAC:TAGCCATG': 1, 'TAGCCATG:CTCTGGAT': 1, 'CTAGCTCA:GCTACTCT': 1, 'CTCTGGAT:TATGGCAC': 1, 'TCTTCGAC:GATCTTGC': 1, 'TCTTCGAC:TACCGGAT': 1, 'TCTTCGAC:GCTACTCT': 1, 'AGGATAGC:TAGCCATG': 1, 'TCGAGAGT:TATGGCAC': 1, 'ATCATGCG:TACCGGAT': 1, 'TACCGGAT:CGATCGAT': 1, 'GATCAAGG:TACCGGAT': 1, 'CGGTAATC:ATCGTGGT': 1, 'ACGATCAG:TCGAGAGT': 1, 'GATCAAGG:AGGATAGC': 1, 'ACGATCAG:TACCGGAT': 1, 'ATCGTGGT:TACCGGAT': 1, 'TACCGGAT:CTAGCTCA': 1, 'CTAGCTCA:AGGATAGC': 1, 'TAGCCATG:TACCGGAT': 1, 'CTAGCTCA:GTCCTAAG': 1, 'TCTTCGAC:TATGGCAC': 1, 'AGGATAGC:ATCATGCG': 1, 'AACAGCGA:CACTTCAC': 1, 'CGGTAATC:TCTTCGAC': 1, 'TCGAGAGT:TGTTCCGT': 1, 'ATCGTGGT:CTCTGGAT': 1, 'TACCGGAT:CGGTAATC': 1, 'TACCGGAT:ACGATCAG': 1, 'CTAGCTCA:GTAGCGTA': 1, 'AGAGTCCA:TCGGATTC': 1, 'GTAGCGTA:TACCGGAT': 1, 'CTCTGGAT:ATCGTGGT': 1, 'TACCGGAT:AACAGCGA': 1, 'TCGAGAGT:AGAGTCCA': 1, 'CTCTGGAT:TAGCCATG': 1, 'TGTTCCGT:TCTTCGAC': 1, 'TACCGGAT:TGTTCCGT': 1, 'TCTTCGAC:TCGAGAGT': 1, 'TCTTCGAC:GTCCTAAG': 1, 'ATCATGCG:AACAGCGA': 1, 'GTAGCGTA:TAGCCATG': 1, 'AGAGTCCA:TCGAGAGT': 1, 'AGAGTCCA:ACGATCAG': 1, 'TACCGGAT:ATCATGCG': 1, 'AACAGCGA:TACCGGAT': 1, 'TATGGCAC:TACCGGAT': 1, 'CTCTGGAT:CTAGCTCA': 1, 'TGTTCCGT:AGAGTCCA': 1, 'ACGATCAG:CTAGCTCA': 1, 'CTAGCTCA:TACCGGAT': 1, 'TCTTCGAC:GATCAAGG': 1, 'GATCAAGG:CTCTGGAT': 1, 'TACCGGAT:CACTTCAC': 1, 'ATCGTGGT:AGGATAGC': 1, 'TCGAGAGT:CTAGCTCA': 1, 'TGTTCCGT:CTCTGGAT': 1, 'GCTACTCT:TGTTCCGT': 1, 'TAGCCATG:TATGGCAC': 1, 'CTAGCTCA:CACTTCAC': 1, 'ATCGTGGT:GATCAAGG': 1, 'CTAGCTCA:TCTTCGAC': 1, 'CTCTGGAT:TCGGATTC': 1, 'TACCGGAT:TAGCCATG': 1, 'GTAGCGTA:TCGAGAGT': 1, 'CGGTAATC:TAGCCATG': 1, 'TACCGGAT:TATGGCAC': 1, 'AGAGTCCA:TATGGCAC': 1, 'TAGCCATG:TCGACAAG': 1, 'CTAGCTCA:ATCATGCG': 1, 'TCGAGAGT:TCTTCGAC': 1, 'TGTTCCGT:TCGACAAG': 1, 'CTCTGGAT:TCGAGAGT': 1, 'GATCTTGC:AACAGCGA': 1, 'TCGAGAGT:CTCTGGAT': 1, 'TACCGGAT:GATCAAGG': 1, 'TCGAGAGT:TACCGGAT': 1, 'TATGGCAC:GTCCTAAG': 1, 'CTCTGGAT:TGTTCCGT': 1, 'ACGATCAG:TATGGCAC': 1, 'AGGATAGC:GATCTTGC': 1, 'CTAGCTCA:TATGGCAC': 1, 'TCGAGAGT:TAGCCATG': 1, 'CTCTGGAT:GTCCTAAG': 1, 'ATCATGCG:TAGCCATG': 1, 'TATGGCAC:CTCTGGAT': 1, 'ATCATGCG:TGTTCCGT': 1, 'TACCGGAT:AGAGTCCA': 1, 'ACGATCAG:AACAGCGA': 1, 'AGAGTCCA:CTCTGGAT': 1, 'GCTACTCT:TCGAGAGT': 1, 'TGTTCCGT:CGATCGAT': 1, 'TCGGATTC:CACTTCAC': 1, 'ATCGTGGT:AGAGTCCA': 1, 'TCGGATTC:TAGCCATG': 1, 'CTCTGGAT:AGAGTCCA': 1, 'GATCTTGC:CTCTGGAT': 1, 'TCGGATTC:TCGACAAG': 1, 'ATCATGCG:TCGACAAG': 1, 'TATGGCAC:CACTTCAC': 1, 'AGGATAGC:CTAGCTCA': 1, 'CACTTCAC:AGAGTCCA': 1, 'AGGATAGC:ACGATCAG': 1, 'TATGGCAC:TAGCCATG': 1, 'TCGACAAG:TCGGATTC': 1, 'AGAGTCCA:CTAGCTCA': 1, 'TATGGCAC:AGAGTCCA': 1, 'ATCGTGGT:CACTTCAC': 1, 'TCGGATTC:GATCAAGG': 1, 'GATCAAGG:CACTTCAC': 1, 'GCTACTCT:TCGACAAG': 1, 'CTCTGGAT:TCTTCGAC': 1, 'TACCGGAT:GTAGCGTA': 1, 'AGAGTCCA:GCTACTCT': 1, 'CGATCGAT:CTCTGGAT': 1, 'AGGATAGC:GATCAAGG': 1, 'CTAGCTCA:AGAGTCCA': 1, 'TCGACAAG:TCGAGAGT': 1, 'ATCATGCG:AGAGTCCA': 1, 'TCGAGAGT:AACAGCGA': 1, 'CGATCGAT:ATCATGCG': 1, 'AGGATAGC:AACAGCGA': 1, 'GCTACTCT:CTCTGGAT': 1, 'TCGGATTC:AACAGCGA': 1, 'TCTTCGAC:CGATCGAT': 1, 'CGATCGAT:AGGATAGC': 1, 'TCTTCGAC:TCGACAAG': 1, 'CACTTCAC:CTAGCTCA': 1, 'CTCTGGAT:GTAGCGTA': 1, 'TCGGATTC:TACCGGAT': 1, 'TCTTCGAC:AGAGTCCA': 1, 'GTAGCGTA:CACTTCAC': 1, 'AACAGCGA:CGGTAATC': 1, 'TGTTCCGT:CACTTCAC': 1, 'GTAGCGTA:CGATCGAT': 1, 'AACAGCGA:TCGAGAGT': 1, 'TACCGGAT:ATCGTGGT': 1, 'AGGATAGC:TGTTCCGT': 1, 'TCGAGAGT:ATCGTGGT': 1, 'ATCATGCG:CGGTAATC': 1, 'GTAGCGTA:GATCTTGC': 1, 'GTAGCGTA:ACGATCAG': 1, 'GATCTTGC:CTAGCTCA': 1, 'TAGCCATG:CACTTCAC': 1, 'TAGCCATG:TGTTCCGT': 1, 'CACTTCAC:GCTACTCT': 1, 'ATCGTGGT:TCTTCGAC': 1, 'AGAGTCCA:ATCGTGGT': 1, 'TAGCCATG:CGGTAATC': 1, 'TACCGGAT:GTCCTAAG': 1, 'CGATCGAT:CTAGCTCA': 1, 'TGTTCCGT:GCTACTCT': 1, 'AACAGCGA:TCTTCGAC': 1, 'TCGGATTC:GCTACTCT': 1, 'CGGTAATC:GTCCTAAG': 1, 'ATCATGCG:ACGATCAG': 1, 'CACTTCAC:TACCGGAT': 1, 'TCTTCGAC:TGTTCCGT': 1, 'TATGGCAC:GCTACTCT': 1, 'AGAGTCCA:AACAGCGA': 1, 'TAGCCATG:TCTTCGAC': 1, 'GTCCTAAG:AACAGCGA': 1, 'TCTTCGAC:CGGTAATC': 1, 'GATCAAGG:TGTTCCGT': 1, 'TCGACAAG:TCTTCGAC': 1, 'CGATCGAT:TCTTCGAC': 1, 'TGTTCCGT:CTAGCTCA': 1, 'ACGATCAG:ATCATGCG': 1, 'GATCAAGG:TAGCCATG': 1, 'TCGAGAGT:AGGATAGC': 1, 'TAGCCATG:CTAGCTCA': 1, 'CACTTCAC:TCGGATTC': 1, 'GTCCTAAG:ACGATCAG': 1, 'TAGCCATG:ACGATCAG': 1, 'ACGATCAG:GTCCTAAG': 1, 'ATCATGCG:CTAGCTCA': 1, 'TGTTCCGT:ATCATGCG': 1, 'ACGATCAG:CTCTGGAT': 1, 'AGAGTCCA:GTAGCGTA': 1, 'TCTTCGAC:AGGATAGC': 1, 'GATCTTGC:CGATCGAT': 1, 'CTAGCTCA:TCGAGAGT': 1, 'CTCTGGAT:CGATCGAT': 1, 'TCGACAAG:CTAGCTCA': 1, 'CTAGCTCA:ATCGTGGT': 1, 'ACGATCAG:TAGCCATG': 1, 'AGGATAGC:TCGGATTC': 1, 'GATCAAGG:GTAGCGTA': 1, 'CTCTGGAT:CACTTCAC': 1, 'AACAGCGA:CTAGCTCA': 1, 'AACAGCGA:CTCTGGAT': 1, 'AACAGCGA:CGATCGAT': 1, 'CGATCGAT:GATCTTGC': 1, 'GATCAAGG:ATCATGCG': 1, 'GTAGCGTA:CTAGCTCA': 1, 'AACAGCGA:ACGATCAG': 1, 'TCGACAAG:TACCGGAT': 1, 'GATCAAGG:CTAGCTCA': 1, 'TCTTCGAC:CTAGCTCA': 1, 'GTCCTAAG:TAGCCATG': 1, 'GTCCTAAG:TCTTCGAC': 1, 'TGTTCCGT:TCGAGAGT': 1, 'GATCTTGC:CGGTAATC': 1, 'CGGTAATC:GATCAAGG': 1, 'CGGTAATC:TGTTCCGT': 1, 'TATGGCAC:TCGAGAGT': 1, 'AGGATAGC:TATGGCAC': 1, 'TATGGCAC:TCGGATTC': 1, 'AACAGCGA:TAGCCATG': 1, 'GCTACTCT:TCTTCGAC': 1, 'TATGGCAC:ATCGTGGT': 1, 'CACTTCAC:TATGGCAC': 1, 'GTCCTAAG:TCGGATTC': 1, 'TATGGCAC:AGGATAGC': 1, 'TCGAGAGT:GCTACTCT': 1, 'TAGCCATG:GATCAAGG': 1, 'TCTTCGAC:ACGATCAG': 1, 'GTAGCGTA:GCTACTCT': 1, 'GATCTTGC:TCTTCGAC': 1, 'TCGGATTC:TCTTCGAC': 1, 'CTCTGGAT:ATCATGCG': 1, 'CGGTAATC:CTCTGGAT': 1, 'CGGTAATC:TCGAGAGT': 1, 'ATCATGCG:CTCTGGAT': 1, 'ACGATCAG:AGAGTCCA': 1, 'TCGGATTC:TCGAGAGT': 1, 'AGGATAGC:TCGAGAGT': 1, 'TCTTCGAC:ATCATGCG': 1, 'AGAGTCCA:TCTTCGAC': 1, 'TCGAGAGT:GATCAAGG': 1, 'GTCCTAAG:GCTACTCT': 1, 'ACGATCAG:ATCGTGGT': 1, 'CTAGCTCA:ACGATCAG': 1, 'TGTTCCGT:GATCAAGG': 1, 'CGGTAATC:ACGATCAG': 1, 'CGGTAATC:AGGATAGC': 1, 'TCGGATTC:CGGTAATC': 1, 'TCTTCGAC:GTAGCGTA': 1, 'AGGATAGC:ATCGTGGT': 1, 'GATCTTGC:TACCGGAT': 1, 'CTAGCTCA:GATCAAGG': 1, 'TCGACAAG:GATCAAGG': 1, 'TCGGATTC:CGATCGAT': 1, 'CGATCGAT:ACGATCAG': 1, 'AACAGCGA:GATCTTGC': 1, 'AACAGCGA:ATCGTGGT': 1, 'ACGATCAG:CACTTCAC': 1, 'GCTACTCT:GTCCTAAG': 1, 'GTAGCGTA:TGTTCCGT': 1, 'TGTTCCGT:TCGGATTC': 1, 'TCTTCGAC:AACAGCGA': 1, 'ACGATCAG:TCTTCGAC': 1, 'AGGATAGC:GTCCTAAG': 1, 'TATGGCAC:ACGATCAG': 1, 'TGTTCCGT:GTCCTAAG': 1, 'TACCGGAT:TCGGATTC': 1, 'TGTTCCGT:ACGATCAG': 1, 'GTAGCGTA:CTCTGGAT': 1, 'GTCCTAAG:TGTTCCGT': 1, 'CGATCGAT:GCTACTCT': 1, 'GTCCTAAG:ATCATGCG': 1, 'TATGGCAC:GTAGCGTA': 1, 'TCGAGAGT:ACGATCAG': 1, 'AGGATAGC:TCGACAAG': 1, 'ATCATGCG:GCTACTCT': 1, 'AACAGCGA:ATCATGCG': 1, 'ACGATCAG:TCGGATTC': 1, 'CTCTGGAT:GATCTTGC': 1, 'CGGTAATC:TCGGATTC': 1, 'CTCTGGAT:ACGATCAG': 1, 'TCGAGAGT:CGATCGAT': 1, 'TGTTCCGT:ATCGTGGT': 1, 'AACAGCGA:GTCCTAAG': 1, 'GTCCTAAG:ATCGTGGT': 1, 'AGGATAGC:CGATCGAT': 1, 'CTAGCTCA:TAGCCATG': 1, 'ACGATCAG:GATCAAGG': 1, 'CGATCGAT:CACTTCAC': 1, 'GTCCTAAG:TCGAGAGT': 1, 'AGAGTCCA:GTCCTAAG': 1, 'GTCCTAAG:AGGATAGC': 1, 'CTCTGGAT:GATCAAGG': 1, 'GTCCTAAG:CTCTGGAT': 1, 'CGGTAATC:TATGGCAC': 1, 'CACTTCAC:CTCTGGAT': 1, 'TCTTCGAC:TAGCCATG': 1, 'GATCAAGG:ATCGTGGT': 1, 'TCGACAAG:GTCCTAAG': 1, 'ATCGTGGT:GTAGCGTA': 1, 'AGGATAGC:CGGTAATC': 1, 'TCGACAAG:CGATCGAT': 1, 'TCGAGAGT:TCGGATTC': 1, 'GCTACTCT:AGGATAGC': 1, 'GATCTTGC:TATGGCAC': 1, 'TGTTCCGT:GTAGCGTA': 1, 'TCGACAAG:TAGCCATG': 1, 'CTAGCTCA:TGTTCCGT': 1, 'TACCGGAT:GATCTTGC': 1, 'GATCTTGC:ACGATCAG': 1, 'TATGGCAC:GATCTTGC': 1, 'AGAGTCCA:ATCATGCG': 1, 'ATCGTGGT:CTAGCTCA': 1, 'ATCGTGGT:TCGAGAGT': 1, 'TCGGATTC:GTCCTAAG': 1, 'TATGGCAC:CTAGCTCA': 1, 'ATCGTGGT:TGTTCCGT': 1, 'ATCATGCG:ATCGTGGT': 1, 'CGGTAATC:CTAGCTCA': 1, 'GCTACTCT:GATCAAGG': 1, 'TCGAGAGT:GTAGCGTA': 1, 'TAGCCATG:GCTACTCT': 1, 'GATCAAGG:CGATCGAT': 1, 'GTAGCGTA:AGAGTCCA': 1, 'TCGGATTC:GATCTTGC': 1, 'CGATCGAT:AACAGCGA': 1, 'CTCTGGAT:CGGTAATC': 1, 'TCGACAAG:TGTTCCGT': 1, 'ACGATCAG:CGGTAATC': 1, 'GTAGCGTA:GATCAAGG': 1, 'CGATCGAT:TATGGCAC': 1, 'GCTACTCT:ACGATCAG': 1, 'GTCCTAAG:GATCAAGG': 1, 'GCTACTCT:GATCTTGC': 1, 'CGGTAATC:TCGACAAG': 1, 'TAGCCATG:ATCATGCG': 1, 'ATCGTGGT:TATGGCAC': 1, 'ATCATGCG:AGGATAGC': 1, 'TGTTCCGT:AGGATAGC': 1, 'TAGCCATG:GATCTTGC': 1, 'AGAGTCCA:TCGACAAG': 1, 'CTCTGGAT:TCGACAAG': 1, 'AACAGCGA:TGTTCCGT': 1, 'GTAGCGTA:TCGGATTC': 1, 'CACTTCAC:AGGATAGC': 1, 'TATGGCAC:ATCATGCG': 1, 'GCTACTCT:TATGGCAC': 1, 'AGGATAGC:AGAGTCCA': 1, 'TATGGCAC:GATCAAGG': 1, 'ATCGTGGT:TAGCCATG': 1, 'TGTTCCGT:CGGTAATC': 1, 'ATCATGCG:TCGAGAGT': 1, 'AGAGTCCA:AGGATAGC': 1, 'ACGATCAG:CGATCGAT': 1, 'ATCATGCG:TCGGATTC': 1, 'CACTTCAC:CGATCGAT': 1, 'AGAGTCCA:CGGTAATC': 1, 'CACTTCAC:ATCATGCG': 1, 'TCGACAAG:AACAGCGA': 1, 'ATCGTGGT:ACGATCAG': 1, 'GTCCTAAG:CGATCGAT': 1, 'TGTTCCGT:GATCTTGC': 1, 'CTAGCTCA:GATCTTGC': 1, 'GATCAAGG:AGAGTCCA': 1, 'CACTTCAC:CGGTAATC': 1, 'TGTTCCGT:AACAGCGA': 1, 'CTAGCTCA:AACAGCGA': 1, 'TCGAGAGT:GTCCTAAG': 1, 'GTAGCGTA:TCTTCGAC': 1, 'GTAGCGTA:ATCGTGGT': 1, 'TCGGATTC:AGGATAGC': 1, 'ATCGTGGT:ATCATGCG': 1, 'TAGCCATG:CGATCGAT': 1, 'GTAGCGTA:AACAGCGA': 1, 'AACAGCGA:GCTACTCT': 1, 'CGGTAATC:AACAGCGA': 1, 'TCGAGAGT:TCGACAAG': 1, 'GCTACTCT:AGAGTCCA': 1, 'AGAGTCCA:TGTTCCGT': 1, 'ATCATGCG:GATCAAGG': 1, 'GTAGCGTA:ATCATGCG': 1, 'TCGACAAG:CTCTGGAT': 1, 'ATCATGCG:GTAGCGTA': 1, 'CGGTAATC:GATCTTGC': 1, 'GTAGCGTA:TATGGCAC': 1, 'TCGACAAG:AGGATAGC': 1, 'CACTTCAC:TGTTCCGT': 1, 'CGGTAATC:CGATCGAT': 1, 'TAGCCATG:TCGGATTC': 1, 'CACTTCAC:GATCTTGC': 1, 'TAGCCATG:AGGATAGC': 1, 'GTCCTAAG:TCGACAAG': 1, 'AACAGCGA:TCGGATTC': 1, 'ATCGTGGT:GTCCTAAG': 1, 'TCGAGAGT:GATCTTGC': 1, 'AGAGTCCA:GATCAAGG': 1, 'ACGATCAG:TGTTCCGT': 1, 'ATCGTGGT:GATCTTGC': 1, 'TAGCCATG:TCGAGAGT': 1, 'GTAGCGTA:CGGTAATC': 1, 'GCTACTCT:TAGCCATG': 1, 'ATCGTGGT:TCGGATTC': 1, 'CGGTAATC:ATCATGCG': 1, 'ATCGTGGT:AACAGCGA': 1, 'GTCCTAAG:GTAGCGTA': 1, 'GCTACTCT:ATCATGCG': 1, 'TCGACAAG:TATGGCAC': 1, 'AACAGCGA:GATCAAGG': 1, 'AACAGCGA:TATGGCAC': 1, 'ATCATGCG:GTCCTAAG': 1, 'TCGAGAGT:CGGTAATC': 1, 'GATCTTGC:GATCAAGG': 1, 'AGAGTCCA:TAGCCATG': 1, 'AACAGCGA:AGAGTCCA': 1, 'GTAGCGTA:GTCCTAAG': 1, 'TAGCCATG:AGAGTCCA': 1, 'TCGACAAG:GATCTTGC': 1, 'GATCTTGC:TCGAGAGT': 1, 'CTAGCTCA:CGATCGAT': 1, 'CGATCGAT:TCGACAAG': 1, 'AACAGCGA:AGGATAGC': 1, 'TAGCCATG:GTAGCGTA': 1, 'GATCTTGC:ATCATGCG': 1, 'ACGATCAG:TCGACAAG': 1, 'CACTTCAC:GTCCTAAG': 1, 'TAGCCATG:GTCCTAAG': 1, 'GATCTTGC:GCTACTCT': 1, 'CTAGCTCA:TCGGATTC': 1, 'CGATCGAT:AGAGTCCA': 1, 'TCGACAAG:ACGATCAG': 1, 'CACTTCAC:TCGAGAGT': 1, 'TCTTCGAC:CACTTCAC': 1, 'TCGACAAG:GTAGCGTA': 1, 'CGATCGAT:TCGAGAGT': 1, 'GCTACTCT:CTAGCTCA': 1, 'TAGCCATG:AACAGCGA': 1, 'TATGGCAC:CGATCGAT': 1, 'AGGATAGC:CACTTCAC': 1, 'TCGGATTC:TATGGCAC': 1, 'ACGATCAG:GATCTTGC': 1, 'CGGTAATC:GTAGCGTA': 1, 'CGATCGAT:GATCAAGG': 1, 'ACGATCAG:AGGATAGC': 1, 'GTAGCGTA:AGGATAGC': 1, 'CGGTAATC:CACTTCAC': 1, 'CGATCGAT:TAGCCATG': 1, 'ATCGTGGT:GCTACTCT': 1, 'TCGAGAGT:ATCATGCG': 1, 'GATCAAGG:GATCTTGC': 1, 'GTCCTAAG:CGGTAATC': 1, 'GCTACTCT:CGGTAATC': 1, 'TCGGATTC:TGTTCCGT': 1, 'CGATCGAT:CGGTAATC': 1, 'GATCAAGG:AACAGCGA': 1, 'GATCAAGG:TCGACAAG': 1, 'TCGGATTC:CTCTGGAT': 1, 'ATCGTGGT:CGGTAATC': 1, 'CGGTAATC:AGAGTCCA': 1, 'GTAGCGTA:TCGACAAG': 1, 'ATCATGCG:GATCTTGC': 1, 'CACTTCAC:AACAGCGA': 1, 'GATCTTGC:TGTTCCGT': 1, 'GATCAAGG:TATGGCAC': 1, 'TATGGCAC:TCGACAAG': 1, 'ATCGTGGT:TCGACAAG': 1, 'GATCTTGC:ATCGTGGT': 1, 'TAGCCATG:ATCGTGGT': 1, 'TCGACAAG:GCTACTCT': 1, 'GTCCTAAG:GATCTTGC': 1, 'GCTACTCT:CGATCGAT': 1, 'TCGGATTC:CTAGCTCA': 1, 'AACAGCGA:GTAGCGTA': 1, 'CGATCGAT:GTCCTAAG': 1, 'ATCGTGGT:CGATCGAT': 1, 'ACGATCAG:GCTACTCT': 1, 'TCGGATTC:ATCGTGGT': 1, 'GTCCTAAG:AGAGTCCA': 1, 'AGAGTCCA:GATCTTGC': 1, 'ATCATGCG:CACTTCAC': 1, 'GATCTTGC:AGGATAGC': 1, 'GATCAAGG:TCGGATTC': 1, 'GATCTTGC:GTCCTAAG': 1, 'GATCAAGG:TCGAGAGT': 1, 'TCGACAAG:ATCGTGGT': 1, 'AGAGTCCA:CACTTCAC': 1, 'ACGATCAG:GTAGCGTA': 1, 'CGATCGAT:ATCGTGGT': 1, 'GATCTTGC:TAGCCATG': 1, 'ATCATGCG:CGATCGAT': 1, 'GCTACTCT:ATCGTGGT': 1, 'GATCTTGC:CACTTCAC': 1, 'AACAGCGA:TCGACAAG': 1, 'TATGGCAC:CGGTAATC': 1, 'GCTACTCT:AACAGCGA': 1, 'GATCTTGC:AGAGTCCA': 1, 'GATCTTGC:TCGACAAG': 1, 'CGGTAATC:GCTACTCT': 1, 'CGATCGAT:TCGGATTC': 1, 'CACTTCAC:GTAGCGTA': 1, 'TCGAGAGT:CACTTCAC': 1, 'GCTACTCT:TCGGATTC': 1, 'CACTTCAC:GATCAAGG': 1, 'CACTTCAC:ATCGTGGT': 1, 'GTCCTAAG:CACTTCAC': 1, 'TCGGATTC:ATCATGCG': 1, 'GATCTTGC:GTAGCGTA': 1, 'GATCTTGC:TCGGATTC': 1, 'TCGACAAG:CGGTAATC': 1, 'TCGGATTC:GTAGCGTA': 1, 'GATCAAGG:CGGTAATC': 1, 'TCGGATTC:AGAGTCCA': 1, 'GCTACTCT:CACTTCAC': 1, 'CACTTCAC:TCGACAAG': 1, 'TCGACAAG:CACTTCAC': 1}
Matched index counts: {'TACCGGAT': 69307073, 'CTCTGGAT': 32163349, 'AGAGTCCA': 10378366, 'GTAGCGTA': 7450201, 'ATCATGCG': 9264615, 'AACAGCGA': 8178191, 'TCGACAAG': 3548541, 'TCGAGAGT': 10658212, 'CGGTAATC': 4498136, 'TAGCCATG': 9852258, 'TCTTCGAC': 39149148, 'CTAGCTCA': 16162895, 'TATGGCAC': 10195805, 'ACGATCAG': 7441721, 'GATCTTGC': 3425453, 'AGGATAGC': 8078057, 'TGTTCCGT': 14786868, 'CACTTCAC': 3833640, 'GCTACTCT': 6610857, 'ATCGTGGT': 6357656, 'CGATCGAT': 5225776, 'GTCCTAAG': 8164223, 'GATCAAGG': 6085915, 'TCGGATTC': 4163314}
Matched index percentages: {'TACCGGAT': 19.079888770369816, 'CTCTGGAT': 8.854408285321545, 'AGAGTCCA': 2.85711198477806, 'GTAGCGTA': 2.0510028810031837, 'ATCATGCG': 2.550501933623712, 'AACAGCGA': 2.2514148681886983, 'TCGACAAG': 0.9768954977668278, 'TCGAGAGT': 2.9341521817119705, 'CGGTAATC': 1.2383142273804608, 'TAGCCATG': 2.7122770972738404, 'TCTTCGAC': 10.777563630406755, 'CTAGCTCA': 4.449563737992029, 'TATGGCAC': 2.8068538592645575, 'ACGATCAG': 2.0486683796345755, 'GATCTTGC': 0.9430099901655, 'AGGATAGC': 2.2238484813910304, 'TGTTCCGT': 4.070750422574342, 'CACTTCAC': 1.055381819192401, 'GCTACTCT': 1.8199356974261585, 'ATCGTGGT': 1.7502307350401924, 'CGATCGAT': 1.4386298613255257, 'GTCCTAAG': 2.2475695480098397, 'GATCAAGG': 1.6754218038601227, 'TCGGATTC': 1.146139414026667}

Realized I didn't have the paranthesis around the whole key in unmatched, so killed and reran. 

Output:
Number of unmatched records: 517612
Number of unknown records: 57748853
Number of total matched records: 304980270
Total Number of Records: 363246735
Unmatched index counts: {'TATGGCAC:TGTTCCGT': 1, 'CTAGCTCA:CGGTAATC': 1, 'CTCTGGAT:AGGATAGC': 1, 'TCGACAAG:ATCATGCG': 1, 'TGTTCCGT:TATGGCAC': 1, 'TACCGGAT:TCGAGAGT': 1, 'GATCAAGG:TCTTCGAC': 1, 'ATCATGCG:TATGGCAC': 1, 'GTCCTAAG:TATGGCAC': 1, 'TCTTCGAC:TCGGATTC': 1, 'TACCGGAT:AGGATAGC': 1, 'GATCAAGG:GCTACTCT': 1, 'CACTTCAC:TCTTCGAC': 1, 'TACCGGAT:TCGACAAG': 1, 'CTCTGGAT:TACCGGAT': 1, 'AGAGTCCA:TACCGGAT': 1, 'CACTTCAC:ACGATCAG': 1, 'TATGGCAC:AACAGCGA': 1, 'CGATCGAT:GTAGCGTA': 1, 'TCTTCGAC:CTCTGGAT': 1, 'CTAGCTCA:CTCTGGAT': 1, 'CGATCGAT:TGTTCCGT': 1, 'AGGATAGC:TCTTCGAC': 1, 'TGTTCCGT:TACCGGAT': 1, 'CTAGCTCA:TCGACAAG': 1, 'TACCGGAT:CTCTGGAT': 1, 'TATGGCAC:TCTTCGAC': 1, 'TCGGATTC:ACGATCAG': 1, 'TACCGGAT:TCTTCGAC': 1, 'CGGTAATC:TACCGGAT': 1, 'GATCAAGG:ACGATCAG': 1, 'CGATCGAT:TACCGGAT': 1, 'TCTTCGAC:ATCGTGGT': 1, 'GTCCTAAG:CTAGCTCA': 1, 'GCTACTCT:GTAGCGTA': 1, 'GCTACTCT:TACCGGAT': 1, 'TGTTCCGT:TAGCCATG': 1, 'AGAGTCCA:CGATCGAT': 1, 'TACCGGAT:GCTACTCT': 1, 'TCGACAAG:AGAGTCCA': 1, 'AGGATAGC:CTCTGGAT': 1, 'GATCAAGG:GTCCTAAG': 1, 'ATCATGCG:TCTTCGAC': 1, 'AGGATAGC:TACCGGAT': 1, 'GTCCTAAG:TACCGGAT': 1, 'AGGATAGC:GTAGCGTA': 1, 'AGGATAGC:GCTACTCT': 1, 'CTCTGGAT:AACAGCGA': 1, 'CTCTGGAT:GCTACTCT': 1, 'CACTTCAC:TAGCCATG': 1, 'TAGCCATG:CTCTGGAT': 1, 'CTAGCTCA:GCTACTCT': 1, 'CTCTGGAT:TATGGCAC': 1, 'TCTTCGAC:GATCTTGC': 1, 'TCTTCGAC:TACCGGAT': 1, 'TCTTCGAC:GCTACTCT': 1, 'AGGATAGC:TAGCCATG': 1, 'TCGAGAGT:TATGGCAC': 1, 'ATCATGCG:TACCGGAT': 1, 'TACCGGAT:CGATCGAT': 1, 'GATCAAGG:TACCGGAT': 1, 'CGGTAATC:ATCGTGGT': 1, 'ACGATCAG:TCGAGAGT': 1, 'GATCAAGG:AGGATAGC': 1, 'ACGATCAG:TACCGGAT': 1, 'ATCGTGGT:TACCGGAT': 1, 'TACCGGAT:CTAGCTCA': 1, 'CTAGCTCA:AGGATAGC': 1, 'TAGCCATG:TACCGGAT': 1, 'CTAGCTCA:GTCCTAAG': 1, 'TCTTCGAC:TATGGCAC': 1, 'AGGATAGC:ATCATGCG': 1, 'AACAGCGA:CACTTCAC': 1, 'CGGTAATC:TCTTCGAC': 1, 'TCGAGAGT:TGTTCCGT': 1, 'ATCGTGGT:CTCTGGAT': 1, 'TACCGGAT:CGGTAATC': 1, 'TACCGGAT:ACGATCAG': 1, 'CTAGCTCA:GTAGCGTA': 1, 'AGAGTCCA:TCGGATTC': 1, 'GTAGCGTA:TACCGGAT': 1, 'CTCTGGAT:ATCGTGGT': 1, 'TACCGGAT:AACAGCGA': 1, 'TCGAGAGT:AGAGTCCA': 1, 'CTCTGGAT:TAGCCATG': 1, 'TGTTCCGT:TCTTCGAC': 1, 'TACCGGAT:TGTTCCGT': 1, 'TCTTCGAC:TCGAGAGT': 1, 'TCTTCGAC:GTCCTAAG': 1, 'ATCATGCG:AACAGCGA': 1, 'GTAGCGTA:TAGCCATG': 1, 'AGAGTCCA:TCGAGAGT': 1, 'AGAGTCCA:ACGATCAG': 1, 'TACCGGAT:ATCATGCG': 1, 'AACAGCGA:TACCGGAT': 1, 'TATGGCAC:TACCGGAT': 1, 'CTCTGGAT:CTAGCTCA': 1, 'TGTTCCGT:AGAGTCCA': 1, 'ACGATCAG:CTAGCTCA': 1, 'CTAGCTCA:TACCGGAT': 1, 'TCTTCGAC:GATCAAGG': 1, 'GATCAAGG:CTCTGGAT': 1, 'TACCGGAT:CACTTCAC': 1, 'ATCGTGGT:AGGATAGC': 1, 'TCGAGAGT:CTAGCTCA': 1, 'TGTTCCGT:CTCTGGAT': 1, 'GCTACTCT:TGTTCCGT': 1, 'TAGCCATG:TATGGCAC': 1, 'CTAGCTCA:CACTTCAC': 1, 'ATCGTGGT:GATCAAGG': 1, 'CTAGCTCA:TCTTCGAC': 1, 'CTCTGGAT:TCGGATTC': 1, 'TACCGGAT:TAGCCATG': 1, 'GTAGCGTA:TCGAGAGT': 1, 'CGGTAATC:TAGCCATG': 1, 'TACCGGAT:TATGGCAC': 1, 'AGAGTCCA:TATGGCAC': 1, 'TAGCCATG:TCGACAAG': 1, 'CTAGCTCA:ATCATGCG': 1, 'TCGAGAGT:TCTTCGAC': 1, 'TGTTCCGT:TCGACAAG': 1, 'CTCTGGAT:TCGAGAGT': 1, 'GATCTTGC:AACAGCGA': 1, 'TCGAGAGT:CTCTGGAT': 1, 'TACCGGAT:GATCAAGG': 1, 'TCGAGAGT:TACCGGAT': 1, 'TATGGCAC:GTCCTAAG': 1, 'CTCTGGAT:TGTTCCGT': 1, 'ACGATCAG:TATGGCAC': 1, 'AGGATAGC:GATCTTGC': 1, 'CTAGCTCA:TATGGCAC': 1, 'TCGAGAGT:TAGCCATG': 1, 'CTCTGGAT:GTCCTAAG': 1, 'ATCATGCG:TAGCCATG': 1, 'TATGGCAC:CTCTGGAT': 1, 'ATCATGCG:TGTTCCGT': 1, 'TACCGGAT:AGAGTCCA': 1, 'ACGATCAG:AACAGCGA': 1, 'AGAGTCCA:CTCTGGAT': 1, 'GCTACTCT:TCGAGAGT': 1, 'TGTTCCGT:CGATCGAT': 1, 'TCGGATTC:CACTTCAC': 1, 'ATCGTGGT:AGAGTCCA': 1, 'TCGGATTC:TAGCCATG': 1, 'CTCTGGAT:AGAGTCCA': 1, 'GATCTTGC:CTCTGGAT': 1, 'TCGGATTC:TCGACAAG': 1, 'ATCATGCG:TCGACAAG': 1, 'TATGGCAC:CACTTCAC': 1, 'AGGATAGC:CTAGCTCA': 1, 'CACTTCAC:AGAGTCCA': 1, 'AGGATAGC:ACGATCAG': 1, 'TATGGCAC:TAGCCATG': 1, 'TCGACAAG:TCGGATTC': 1, 'AGAGTCCA:CTAGCTCA': 1, 'TATGGCAC:AGAGTCCA': 1, 'ATCGTGGT:CACTTCAC': 1, 'TCGGATTC:GATCAAGG': 1, 'GATCAAGG:CACTTCAC': 1, 'GCTACTCT:TCGACAAG': 1, 'CTCTGGAT:TCTTCGAC': 1, 'TACCGGAT:GTAGCGTA': 1, 'AGAGTCCA:GCTACTCT': 1, 'CGATCGAT:CTCTGGAT': 1, 'AGGATAGC:GATCAAGG': 1, 'CTAGCTCA:AGAGTCCA': 1, 'TCGACAAG:TCGAGAGT': 1, 'ATCATGCG:AGAGTCCA': 1, 'TCGAGAGT:AACAGCGA': 1, 'CGATCGAT:ATCATGCG': 1, 'AGGATAGC:AACAGCGA': 1, 'GCTACTCT:CTCTGGAT': 1, 'TCGGATTC:AACAGCGA': 1, 'TCTTCGAC:CGATCGAT': 1, 'CGATCGAT:AGGATAGC': 1, 'TCTTCGAC:TCGACAAG': 1, 'CACTTCAC:CTAGCTCA': 1, 'CTCTGGAT:GTAGCGTA': 1, 'TCGGATTC:TACCGGAT': 1, 'TCTTCGAC:AGAGTCCA': 1, 'GTAGCGTA:CACTTCAC': 1, 'AACAGCGA:CGGTAATC': 1, 'TGTTCCGT:CACTTCAC': 1, 'GTAGCGTA:CGATCGAT': 1, 'AACAGCGA:TCGAGAGT': 1, 'TACCGGAT:ATCGTGGT': 1, 'AGGATAGC:TGTTCCGT': 1, 'TCGAGAGT:ATCGTGGT': 1, 'ATCATGCG:CGGTAATC': 1, 'GTAGCGTA:GATCTTGC': 1, 'GTAGCGTA:ACGATCAG': 1, 'GATCTTGC:CTAGCTCA': 1, 'TAGCCATG:CACTTCAC': 1, 'TAGCCATG:TGTTCCGT': 1, 'CACTTCAC:GCTACTCT': 1, 'ATCGTGGT:TCTTCGAC': 1, 'AGAGTCCA:ATCGTGGT': 1, 'TAGCCATG:CGGTAATC': 1, 'TACCGGAT:GTCCTAAG': 1, 'CGATCGAT:CTAGCTCA': 1, 'TGTTCCGT:GCTACTCT': 1, 'AACAGCGA:TCTTCGAC': 1, 'TCGGATTC:GCTACTCT': 1, 'CGGTAATC:GTCCTAAG': 1, 'ATCATGCG:ACGATCAG': 1, 'CACTTCAC:TACCGGAT': 1, 'TCTTCGAC:TGTTCCGT': 1, 'TATGGCAC:GCTACTCT': 1, 'AGAGTCCA:AACAGCGA': 1, 'TAGCCATG:TCTTCGAC': 1, 'GTCCTAAG:AACAGCGA': 1, 'TCTTCGAC:CGGTAATC': 1, 'GATCAAGG:TGTTCCGT': 1, 'TCGACAAG:TCTTCGAC': 1, 'CGATCGAT:TCTTCGAC': 1, 'TGTTCCGT:CTAGCTCA': 1, 'ACGATCAG:ATCATGCG': 1, 'GATCAAGG:TAGCCATG': 1, 'TCGAGAGT:AGGATAGC': 1, 'TAGCCATG:CTAGCTCA': 1, 'CACTTCAC:TCGGATTC': 1, 'GTCCTAAG:ACGATCAG': 1, 'TAGCCATG:ACGATCAG': 1, 'ACGATCAG:GTCCTAAG': 1, 'ATCATGCG:CTAGCTCA': 1, 'TGTTCCGT:ATCATGCG': 1, 'ACGATCAG:CTCTGGAT': 1, 'AGAGTCCA:GTAGCGTA': 1, 'TCTTCGAC:AGGATAGC': 1, 'GATCTTGC:CGATCGAT': 1, 'CTAGCTCA:TCGAGAGT': 1, 'CTCTGGAT:CGATCGAT': 1, 'TCGACAAG:CTAGCTCA': 1, 'CTAGCTCA:ATCGTGGT': 1, 'ACGATCAG:TAGCCATG': 1, 'AGGATAGC:TCGGATTC': 1, 'GATCAAGG:GTAGCGTA': 1, 'CTCTGGAT:CACTTCAC': 1, 'AACAGCGA:CTAGCTCA': 1, 'AACAGCGA:CTCTGGAT': 1, 'AACAGCGA:CGATCGAT': 1, 'CGATCGAT:GATCTTGC': 1, 'GATCAAGG:ATCATGCG': 1, 'GTAGCGTA:CTAGCTCA': 1, 'AACAGCGA:ACGATCAG': 1, 'TCGACAAG:TACCGGAT': 1, 'GATCAAGG:CTAGCTCA': 1, 'TCTTCGAC:CTAGCTCA': 1, 'GTCCTAAG:TAGCCATG': 1, 'GTCCTAAG:TCTTCGAC': 1, 'TGTTCCGT:TCGAGAGT': 1, 'GATCTTGC:CGGTAATC': 1, 'CGGTAATC:GATCAAGG': 1, 'CGGTAATC:TGTTCCGT': 1, 'TATGGCAC:TCGAGAGT': 1, 'AGGATAGC:TATGGCAC': 1, 'TATGGCAC:TCGGATTC': 1, 'AACAGCGA:TAGCCATG': 1, 'GCTACTCT:TCTTCGAC': 1, 'TATGGCAC:ATCGTGGT': 1, 'CACTTCAC:TATGGCAC': 1, 'GTCCTAAG:TCGGATTC': 1, 'TATGGCAC:AGGATAGC': 1, 'TCGAGAGT:GCTACTCT': 1, 'TAGCCATG:GATCAAGG': 1, 'TCTTCGAC:ACGATCAG': 1, 'GTAGCGTA:GCTACTCT': 1, 'GATCTTGC:TCTTCGAC': 1, 'TCGGATTC:TCTTCGAC': 1, 'CTCTGGAT:ATCATGCG': 1, 'CGGTAATC:CTCTGGAT': 1, 'CGGTAATC:TCGAGAGT': 1, 'ATCATGCG:CTCTGGAT': 1, 'ACGATCAG:AGAGTCCA': 1, 'TCGGATTC:TCGAGAGT': 1, 'AGGATAGC:TCGAGAGT': 1, 'TCTTCGAC:ATCATGCG': 1, 'AGAGTCCA:TCTTCGAC': 1, 'TCGAGAGT:GATCAAGG': 1, 'GTCCTAAG:GCTACTCT': 1, 'ACGATCAG:ATCGTGGT': 1, 'CTAGCTCA:ACGATCAG': 1, 'TGTTCCGT:GATCAAGG': 1, 'CGGTAATC:ACGATCAG': 1, 'CGGTAATC:AGGATAGC': 1, 'TCGGATTC:CGGTAATC': 1, 'TCTTCGAC:GTAGCGTA': 1, 'AGGATAGC:ATCGTGGT': 1, 'GATCTTGC:TACCGGAT': 1, 'CTAGCTCA:GATCAAGG': 1, 'TCGACAAG:GATCAAGG': 1, 'TCGGATTC:CGATCGAT': 1, 'CGATCGAT:ACGATCAG': 1, 'AACAGCGA:GATCTTGC': 1, 'AACAGCGA:ATCGTGGT': 1, 'ACGATCAG:CACTTCAC': 1, 'GCTACTCT:GTCCTAAG': 1, 'GTAGCGTA:TGTTCCGT': 1, 'TGTTCCGT:TCGGATTC': 1, 'TCTTCGAC:AACAGCGA': 1, 'ACGATCAG:TCTTCGAC': 1, 'AGGATAGC:GTCCTAAG': 1, 'TATGGCAC:ACGATCAG': 1, 'TGTTCCGT:GTCCTAAG': 1, 'TACCGGAT:TCGGATTC': 1, 'TGTTCCGT:ACGATCAG': 1, 'GTAGCGTA:CTCTGGAT': 1, 'GTCCTAAG:TGTTCCGT': 1, 'CGATCGAT:GCTACTCT': 1, 'GTCCTAAG:ATCATGCG': 1, 'TATGGCAC:GTAGCGTA': 1, 'TCGAGAGT:ACGATCAG': 1, 'AGGATAGC:TCGACAAG': 1, 'ATCATGCG:GCTACTCT': 1, 'AACAGCGA:ATCATGCG': 1, 'ACGATCAG:TCGGATTC': 1, 'CTCTGGAT:GATCTTGC': 1, 'CGGTAATC:TCGGATTC': 1, 'CTCTGGAT:ACGATCAG': 1, 'TCGAGAGT:CGATCGAT': 1, 'TGTTCCGT:ATCGTGGT': 1, 'AACAGCGA:GTCCTAAG': 1, 'GTCCTAAG:ATCGTGGT': 1, 'AGGATAGC:CGATCGAT': 1, 'CTAGCTCA:TAGCCATG': 1, 'ACGATCAG:GATCAAGG': 1, 'CGATCGAT:CACTTCAC': 1, 'GTCCTAAG:TCGAGAGT': 1, 'AGAGTCCA:GTCCTAAG': 1, 'GTCCTAAG:AGGATAGC': 1, 'CTCTGGAT:GATCAAGG': 1, 'GTCCTAAG:CTCTGGAT': 1, 'CGGTAATC:TATGGCAC': 1, 'CACTTCAC:CTCTGGAT': 1, 'TCTTCGAC:TAGCCATG': 1, 'GATCAAGG:ATCGTGGT': 1, 'TCGACAAG:GTCCTAAG': 1, 'ATCGTGGT:GTAGCGTA': 1, 'AGGATAGC:CGGTAATC': 1, 'TCGACAAG:CGATCGAT': 1, 'TCGAGAGT:TCGGATTC': 1, 'GCTACTCT:AGGATAGC': 1, 'GATCTTGC:TATGGCAC': 1, 'TGTTCCGT:GTAGCGTA': 1, 'TCGACAAG:TAGCCATG': 1, 'CTAGCTCA:TGTTCCGT': 1, 'TACCGGAT:GATCTTGC': 1, 'GATCTTGC:ACGATCAG': 1, 'TATGGCAC:GATCTTGC': 1, 'AGAGTCCA:ATCATGCG': 1, 'ATCGTGGT:CTAGCTCA': 1, 'ATCGTGGT:TCGAGAGT': 1, 'TCGGATTC:GTCCTAAG': 1, 'TATGGCAC:CTAGCTCA': 1, 'ATCGTGGT:TGTTCCGT': 1, 'ATCATGCG:ATCGTGGT': 1, 'CGGTAATC:CTAGCTCA': 1, 'GCTACTCT:GATCAAGG': 1, 'TCGAGAGT:GTAGCGTA': 1, 'TAGCCATG:GCTACTCT': 1, 'GATCAAGG:CGATCGAT': 1, 'GTAGCGTA:AGAGTCCA': 1, 'TCGGATTC:GATCTTGC': 1, 'CGATCGAT:AACAGCGA': 1, 'CTCTGGAT:CGGTAATC': 1, 'TCGACAAG:TGTTCCGT': 1, 'ACGATCAG:CGGTAATC': 1, 'GTAGCGTA:GATCAAGG': 1, 'CGATCGAT:TATGGCAC': 1, 'GCTACTCT:ACGATCAG': 1, 'GTCCTAAG:GATCAAGG': 1, 'GCTACTCT:GATCTTGC': 1, 'CGGTAATC:TCGACAAG': 1, 'TAGCCATG:ATCATGCG': 1, 'ATCGTGGT:TATGGCAC': 1, 'ATCATGCG:AGGATAGC': 1, 'TGTTCCGT:AGGATAGC': 1, 'TAGCCATG:GATCTTGC': 1, 'AGAGTCCA:TCGACAAG': 1, 'CTCTGGAT:TCGACAAG': 1, 'AACAGCGA:TGTTCCGT': 1, 'GTAGCGTA:TCGGATTC': 1, 'CACTTCAC:AGGATAGC': 1, 'TATGGCAC:ATCATGCG': 1, 'GCTACTCT:TATGGCAC': 1, 'AGGATAGC:AGAGTCCA': 1, 'TATGGCAC:GATCAAGG': 1, 'ATCGTGGT:TAGCCATG': 1, 'TGTTCCGT:CGGTAATC': 1, 'ATCATGCG:TCGAGAGT': 1, 'AGAGTCCA:AGGATAGC': 1, 'ACGATCAG:CGATCGAT': 1, 'ATCATGCG:TCGGATTC': 1, 'CACTTCAC:CGATCGAT': 1, 'AGAGTCCA:CGGTAATC': 1, 'CACTTCAC:ATCATGCG': 1, 'TCGACAAG:AACAGCGA': 1, 'ATCGTGGT:ACGATCAG': 1, 'GTCCTAAG:CGATCGAT': 1, 'TGTTCCGT:GATCTTGC': 1, 'CTAGCTCA:GATCTTGC': 1, 'GATCAAGG:AGAGTCCA': 1, 'CACTTCAC:CGGTAATC': 1, 'TGTTCCGT:AACAGCGA': 1, 'CTAGCTCA:AACAGCGA': 1, 'TCGAGAGT:GTCCTAAG': 1, 'GTAGCGTA:TCTTCGAC': 1, 'GTAGCGTA:ATCGTGGT': 1, 'TCGGATTC:AGGATAGC': 1, 'ATCGTGGT:ATCATGCG': 1, 'TAGCCATG:CGATCGAT': 1, 'GTAGCGTA:AACAGCGA': 1, 'AACAGCGA:GCTACTCT': 1, 'CGGTAATC:AACAGCGA': 1, 'TCGAGAGT:TCGACAAG': 1, 'GCTACTCT:AGAGTCCA': 1, 'AGAGTCCA:TGTTCCGT': 1, 'ATCATGCG:GATCAAGG': 1, 'GTAGCGTA:ATCATGCG': 1, 'TCGACAAG:CTCTGGAT': 1, 'ATCATGCG:GTAGCGTA': 1, 'CGGTAATC:GATCTTGC': 1, 'GTAGCGTA:TATGGCAC': 1, 'TCGACAAG:AGGATAGC': 1, 'CACTTCAC:TGTTCCGT': 1, 'CGGTAATC:CGATCGAT': 1, 'TAGCCATG:TCGGATTC': 1, 'CACTTCAC:GATCTTGC': 1, 'TAGCCATG:AGGATAGC': 1, 'GTCCTAAG:TCGACAAG': 1, 'AACAGCGA:TCGGATTC': 1, 'ATCGTGGT:GTCCTAAG': 1, 'TCGAGAGT:GATCTTGC': 1, 'AGAGTCCA:GATCAAGG': 1, 'ACGATCAG:TGTTCCGT': 1, 'ATCGTGGT:GATCTTGC': 1, 'TAGCCATG:TCGAGAGT': 1, 'GTAGCGTA:CGGTAATC': 1, 'GCTACTCT:TAGCCATG': 1, 'ATCGTGGT:TCGGATTC': 1, 'CGGTAATC:ATCATGCG': 1, 'ATCGTGGT:AACAGCGA': 1, 'GTCCTAAG:GTAGCGTA': 1, 'GCTACTCT:ATCATGCG': 1, 'TCGACAAG:TATGGCAC': 1, 'AACAGCGA:GATCAAGG': 1, 'AACAGCGA:TATGGCAC': 1, 'ATCATGCG:GTCCTAAG': 1, 'TCGAGAGT:CGGTAATC': 1, 'GATCTTGC:GATCAAGG': 1, 'AGAGTCCA:TAGCCATG': 1, 'AACAGCGA:AGAGTCCA': 1, 'GTAGCGTA:GTCCTAAG': 1, 'TAGCCATG:AGAGTCCA': 1, 'TCGACAAG:GATCTTGC': 1, 'GATCTTGC:TCGAGAGT': 1, 'CTAGCTCA:CGATCGAT': 1, 'CGATCGAT:TCGACAAG': 1, 'AACAGCGA:AGGATAGC': 1, 'TAGCCATG:GTAGCGTA': 1, 'GATCTTGC:ATCATGCG': 1, 'ACGATCAG:TCGACAAG': 1, 'CACTTCAC:GTCCTAAG': 1, 'TAGCCATG:GTCCTAAG': 1, 'GATCTTGC:GCTACTCT': 1, 'CTAGCTCA:TCGGATTC': 1, 'CGATCGAT:AGAGTCCA': 1, 'TCGACAAG:ACGATCAG': 1, 'CACTTCAC:TCGAGAGT': 1, 'TCTTCGAC:CACTTCAC': 1, 'TCGACAAG:GTAGCGTA': 1, 'CGATCGAT:TCGAGAGT': 1, 'GCTACTCT:CTAGCTCA': 1, 'TAGCCATG:AACAGCGA': 1, 'TATGGCAC:CGATCGAT': 1, 'AGGATAGC:CACTTCAC': 1, 'TCGGATTC:TATGGCAC': 1, 'ACGATCAG:GATCTTGC': 1, 'CGGTAATC:GTAGCGTA': 1, 'CGATCGAT:GATCAAGG': 1, 'ACGATCAG:AGGATAGC': 1, 'GTAGCGTA:AGGATAGC': 1, 'CGGTAATC:CACTTCAC': 1, 'CGATCGAT:TAGCCATG': 1, 'ATCGTGGT:GCTACTCT': 1, 'TCGAGAGT:ATCATGCG': 1, 'GATCAAGG:GATCTTGC': 1, 'GTCCTAAG:CGGTAATC': 1, 'GCTACTCT:CGGTAATC': 1, 'TCGGATTC:TGTTCCGT': 1, 'CGATCGAT:CGGTAATC': 1, 'GATCAAGG:AACAGCGA': 1, 'GATCAAGG:TCGACAAG': 1, 'TCGGATTC:CTCTGGAT': 1, 'ATCGTGGT:CGGTAATC': 1, 'CGGTAATC:AGAGTCCA': 1, 'GTAGCGTA:TCGACAAG': 1, 'ATCATGCG:GATCTTGC': 1, 'CACTTCAC:AACAGCGA': 1, 'GATCTTGC:TGTTCCGT': 1, 'GATCAAGG:TATGGCAC': 1, 'TATGGCAC:TCGACAAG': 1, 'ATCGTGGT:TCGACAAG': 1, 'GATCTTGC:ATCGTGGT': 1, 'TAGCCATG:ATCGTGGT': 1, 'TCGACAAG:GCTACTCT': 1, 'GTCCTAAG:GATCTTGC': 1, 'GCTACTCT:CGATCGAT': 1, 'TCGGATTC:CTAGCTCA': 1, 'AACAGCGA:GTAGCGTA': 1, 'CGATCGAT:GTCCTAAG': 1, 'ATCGTGGT:CGATCGAT': 1, 'ACGATCAG:GCTACTCT': 1, 'TCGGATTC:ATCGTGGT': 1, 'GTCCTAAG:AGAGTCCA': 1, 'AGAGTCCA:GATCTTGC': 1, 'ATCATGCG:CACTTCAC': 1, 'GATCTTGC:AGGATAGC': 1, 'GATCAAGG:TCGGATTC': 1, 'GATCTTGC:GTCCTAAG': 1, 'GATCAAGG:TCGAGAGT': 1, 'TCGACAAG:ATCGTGGT': 1, 'AGAGTCCA:CACTTCAC': 1, 'ACGATCAG:GTAGCGTA': 1, 'CGATCGAT:ATCGTGGT': 1, 'GATCTTGC:TAGCCATG': 1, 'ATCATGCG:CGATCGAT': 1, 'GCTACTCT:ATCGTGGT': 1, 'GATCTTGC:CACTTCAC': 1, 'AACAGCGA:TCGACAAG': 1, 'TATGGCAC:CGGTAATC': 1, 'GCTACTCT:AACAGCGA': 1, 'GATCTTGC:AGAGTCCA': 1, 'GATCTTGC:TCGACAAG': 1, 'CGGTAATC:GCTACTCT': 1, 'CGATCGAT:TCGGATTC': 1, 'CACTTCAC:GTAGCGTA': 1, 'TCGAGAGT:CACTTCAC': 1, 'GCTACTCT:TCGGATTC': 1, 'CACTTCAC:GATCAAGG': 1, 'CACTTCAC:ATCGTGGT': 1, 'GTCCTAAG:CACTTCAC': 1, 'TCGGATTC:ATCATGCG': 1, 'GATCTTGC:GTAGCGTA': 1, 'GATCTTGC:TCGGATTC': 1, 'TCGACAAG:CGGTAATC': 1, 'TCGGATTC:GTAGCGTA': 1, 'GATCAAGG:CGGTAATC': 1, 'TCGGATTC:AGAGTCCA': 1, 'GCTACTCT:CACTTCAC': 1, 'CACTTCAC:TCGACAAG': 1, 'TCGACAAG:CACTTCAC': 1}
Matched index counts: {'TACCGGAT': 69307073, 'CTCTGGAT': 32163349, 'AGAGTCCA': 10378366, 'GTAGCGTA': 7450201, 'ATCATGCG': 9264615, 'AACAGCGA': 8178191, 'TCGACAAG': 3548541, 'TCGAGAGT': 10658212, 'CGGTAATC': 4498136, 'TAGCCATG': 9852258, 'TCTTCGAC': 39149148, 'CTAGCTCA': 16162895, 'TATGGCAC': 10195805, 'ACGATCAG': 7441721, 'GATCTTGC': 3425453, 'AGGATAGC': 8078057, 'TGTTCCGT': 14786868, 'CACTTCAC': 3833640, 'GCTACTCT': 6610857, 'ATCGTGGT': 6357656, 'CGATCGAT': 5225776, 'GTCCTAAG': 8164223, 'GATCAAGG': 6085915, 'TCGGATTC': 4163314}
Matched index percentages: {'TACCGGAT': 19.079888770369816, 'CTCTGGAT': 8.854408285321545, 'AGAGTCCA': 2.85711198477806, 'GTAGCGTA': 2.0510028810031837, 'ATCATGCG': 2.550501933623712, 'AACAGCGA': 2.2514148681886983, 'TCGACAAG': 0.9768954977668278, 'TCGAGAGT': 2.9341521817119705, 'CGGTAATC': 1.2383142273804608, 'TAGCCATG': 2.7122770972738404, 'TCTTCGAC': 10.777563630406755, 'CTAGCTCA': 4.449563737992029, 'TATGGCAC': 2.8068538592645575, 'ACGATCAG': 2.0486683796345755, 'GATCTTGC': 0.9430099901655, 'AGGATAGC': 2.2238484813910304, 'TGTTCCGT': 4.070750422574342, 'CACTTCAC': 1.055381819192401, 'GCTACTCT': 1.8199356974261585, 'ATCGTGGT': 1.7502307350401924, 'CGATCGAT': 1.4386298613255257, 'GTCCTAAG': 2.2475695480098397, 'GATCAAGG': 1.6754218038601227, 'TCGGATTC': 1.146139414026667}

Unmatched dictionary is still wrong. I'm just going to count if index 1 is in the dictionary instead.
changed:
#count how many mismatches
if str(index1_record[1] + rev_comp_index2) in unmatched_index_dict:
    unmatched_index_dict[str(index1_record[1] + ':' + rev_comp_index2)] += 1
else:
    unmatched_index_dict[str(index1_record[1] + ':' + rev_comp_index2)] = 1
to:
#count how many mismatches
if index1_record[1] in unmatched_index_dict:
    unmatched_index_dict[index1_record[1]] += 1
else:
    unmatched_index_dict[index1_record[1]] = 1

Reran on talapas.
Command being timed: "./demux_alg.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -index /projects/bgmp/shared/2017_sequencing/indexes.txt -outdir /projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results"
	User time (seconds): 12403.66
	System time (seconds): 73.12
	Percent of CPU this job got: 96%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:36:18
    
Output:
Number of unmatched records: 517612
Number of unknown records: 57748853
Number of total matched records: 304980270
Total Number of Records: 363246735
Unmatched index counts: {'TATGGCAC': 97103, 'CTAGCTCA': 31880, 'CTCTGGAT': 34554, 'TCGACAAG': 10588, 'TGTTCCGT': 92569, 'TACCGGAT': 55532, 'GATCAAGG': 20717, 'ATCATGCG': 10456, 'GTCCTAAG': 14018, 'TCTTCGAC': 38630, 'CACTTCAC': 7172, 'AGAGTCCA': 9889, 'CGATCGAT': 6368, 'AGGATAGC': 8901, 'TCGGATTC': 5516, 'CGGTAATC': 9848, 'GCTACTCT': 6990, 'TAGCCATG': 7885, 'TCGAGAGT': 11410, 'ACGATCAG': 8142, 'ATCGTGGT': 7218, 'AACAGCGA': 9747, 'GTAGCGTA': 8788, 'GATCTTGC': 3691}
Matched index counts: {'TACCGGAT': 69307073, 'CTCTGGAT': 32163349, 'AGAGTCCA': 10378366, 'GTAGCGTA': 7450201, 'ATCATGCG': 9264615, 'AACAGCGA': 8178191, 'TCGACAAG': 3548541, 'TCGAGAGT': 10658212, 'CGGTAATC': 4498136, 'TAGCCATG': 9852258, 'TCTTCGAC': 39149148, 'CTAGCTCA': 16162895, 'TATGGCAC': 10195805, 'ACGATCAG': 7441721, 'GATCTTGC': 3425453, 'AGGATAGC': 8078057, 'TGTTCCGT': 14786868, 'CACTTCAC': 3833640, 'GCTACTCT': 6610857, 'ATCGTGGT': 6357656, 'CGATCGAT': 5225776, 'GTCCTAAG': 8164223, 'GATCAAGG': 6085915, 'TCGGATTC': 4163314}
Matched index percentages: {'TACCGGAT': 19.079888770369816, 'CTCTGGAT': 8.854408285321545, 'AGAGTCCA': 2.85711198477806, 'GTAGCGTA': 2.0510028810031837, 'ATCATGCG': 2.550501933623712, 'AACAGCGA': 2.2514148681886983, 'TCGACAAG': 0.9768954977668278, 'TCGAGAGT': 2.9341521817119705, 'CGGTAATC': 1.2383142273804608, 'TAGCCATG': 2.7122770972738404, 'TCTTCGAC': 10.777563630406755, 'CTAGCTCA': 4.449563737992029, 'TATGGCAC': 2.8068538592645575, 'ACGATCAG': 2.0486683796345755, 'GATCTTGC': 0.9430099901655, 'AGGATAGC': 2.2238484813910304, 'TGTTCCGT': 4.070750422574342, 'CACTTCAC': 1.055381819192401, 'GCTACTCT': 1.8199356974261585, 'ATCGTGGT': 1.7502307350401924, 'CGATCGAT': 1.4386298613255257, 'GTCCTAAG': 2.2475695480098397, 'GATCAAGG': 1.6754218038601227, 'TCGGATTC': 1.146139414026667}
