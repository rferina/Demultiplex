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
```zcat 1294_S1_L008_R2_001.fastq.gz | sed -n 4~4p | grep -E "[a-z]+"```
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
# def populate_array(file): 
#     """Opens a FASTQ file and decodes Phred quality scores to numbers
#     accounting for Phred+33. Sums the quality scores for each position
#     and counts the total number of lines in the FASTQ file.
#     Returns the array and line count."""
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
```zcat 1294_S1_L008_R2_001.fastq.gz | sed -n 2~4p | grep "N" | wc -l```
```
3976613
```zcat 1294_S1_L008_R3_001.fastq.gz | sed -n 2~4p | grep "N" | wc -l
```
3328051

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

