Rachel Ferina 
Lab Notebook: Demux Part 1

Python version: Python 3.10.4
Environment: 

**26 July 2022**

Created public Demultiplex repository by using the template from https://github.com/Leslie-C/Demultiplexing. 
Cloned Demultiplex repository onto talapas in the directory demux that I created, /projects/bgmp/rferina/bioinfo/Bi622/demux

The data is located on talapas in /projects/bgmp/shared/2017_sequencing/
It should not be unzipped, as the files are very large.
    1294_S1_L008_R1_001.fastq.gz
    1294_S1_L008_R2_001.fastq.gz
    1294_S1_L008_R3_001.fastq.gz
    1294_S1_L008_R4_001.fastq.gz

Worked on Part 1:

```command```
| File name | Contents | 
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 | 
| 1294_S1_L008_R4_001.fastq.gz | Read 2 |


Read length head -2 | tail -1

If the data was Phred+64 encoded, the quality score lines would include lowercase letters that translate to higher Phred scores? Ran on the index files to reduce runtime as they have shorter sequences.
```zcat 1294_S1_L008_R2_001.fastq.gz | sed -n 4~4p | grep -E "[a-z]+"```
```
No output was found, so this data is Phred+33 encoded.
Also has # 

confirmed can do numpy on a list
# lst = [4, 5, 7, 8, 9]
# # lst = numpy.ones(5)
# print(numpy.mean(lst))

sed to seq line, grep N 

decided not to use init list, as it is only needed for jupyter notebook

Made a test file with 3 records to see if the histograms were being produced correctly. Had to gzip it to be consistent with the real files.
    gzip test.fq
Now its called test.fq.gz

chmod 755 processing.srun 
forgot ./ in sbatch script

N: amount of indexes that have N's in index files
```zcat 1294_S1_L008_R2_001.fastq.gz | sed -n 2~4p | grep "N" | wc -l```
```
3976613
```zcat 1294_S1_L008_R3_001.fastq.gz | sed -n 2~4p | grep "N" | wc -l
```
3328051
rev comp flip string, translate

graph r3=2 is supposed to be r2

quality scores below 30 excluded.

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


**27 July 2022**

Worked on part 2.

Diagrammed pseudocode options on whiteboard with Lisa, Justine, Jessica, and Kaitlyn.

decided against
dict with keys of forward index, vals rev comp

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


for header line
fh.write(f"{read1[0]} {index1[1]} {index2}\n")
fh.write(f"{read1[1]}\n{read1[2]}\n{read1[3]}\n") # lines 2-4
or str: index1 or index2[1]

first 8 records


decided on set rather than dict

**28 July 2022**

wrote unit tests
index 1 is forward, index 2 is rev comp

changed indexes, some not pass, some not match

can use qual score function if decide want to do average of qscores rather than idividual for low qual

test unknown , unknown N, unknown unmatched

make index1, right index2