#!/bin/sh
#$ -cwd
#$ -R y
#$ -l mf=1G,h_vmem=2G,h_fsize=3G
#$ -m e -M leslie.myint@gmail.com
module load bowtie
module load samtools
cd ~/statgenomics/lab_seqalign
bowtie --sam ~/bowtie_indexes/s_cerevisiae example.fastq | samtools view -b - > example_batchjob.bam