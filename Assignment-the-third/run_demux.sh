#!/bin/bash

#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=demux_%j      ### Job Name
#SBATCH --output=demux_%j.out         ### File in which to store job output
#SBATCH --error=demux_%j.err          ### File in which to store job error messages
#SBATCH --time=0-10:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --account=bgmp      ### Account used for job submission

results='/projects/bgmp/rferina/bioinfo/Bi622/demux/Demultiplex/Assignment-the-third/results'
data_dir='/projects/bgmp/shared/2017_sequencing'
read1='1294_S1_L008_R1_001.fastq.gz'
index1='1294_S1_L008_R2_001.fastq.gz'
index2='1294_S1_L008_R3_001.fastq.gz'
read2='1294_S1_L008_R4_001.fastq.gz'
indexes='indexes.txt'

/usr/bin/time -v ./demux_alg.py -r1 $data_dir/$read1 -i1 $data_dir/$index1 -i2 $data_dir/$index2 -r2 $data_dir/$read2 -index $data_dir/$indexes -outdir $results 

/usr/bin/time -v gzip $results/*.fq

