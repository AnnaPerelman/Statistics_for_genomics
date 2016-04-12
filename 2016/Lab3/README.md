Cheat sheet for Unix: [http://www.rain.org/~mkummel/unix.html](http://www.rain.org/~mkummel/unix.html)

Tutorials for Unix: [http://www.ee.surrey.ac.uk/Teaching/Unix/](http://www.ee.surrey.ac.uk/Teaching/Unix/)

### Getting started

To login to the cluster:

	ssh lmyint1@jhpce01.jhsph.edu

Basic commands:

	mkdir <nameOfYourDirectory> # To create a directory 
	cd <nameOfYourDirectory> # To access a directory
	touch <newFile> # To create a new file

To see what software modules are available (hit enter or space to see more, q to exit the menu):

	module avail

We first get on a node using the qrsh command. By specifying rnet, you login to a node optimized for downloading files.

	qrsh -l rnet mem_free=1G

### Some downloads:

Let's first set up our folders.
	
	mkdir statgenomics
	cd statgenomics
	mkdir lab_seqalign
	cd lab_seqalign

To download the SRA experiment in the lab_seqalign folder:

	curl -O ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP%2FSRP014%2FSRP014134/SRR518875/SRR518875.sra

Log off of this node and log in to a standard one:

	exit
	qrsh -l mem_free=1G

### Converting the SRA experiment to fastq file:

Load SRA toolkit:

	module load sratoolkit
	
We are ready to extract the fastq file from the SRA experiment:

	cd statgenomics/lab_seqalign
	fastq-dump SRR518875.sra

We now have the file SRR518875.fastq in the lab_seqalign directory. 

	head -n 8 SRR518875.fastq
	wc -l SRR518875.fastq # Number of reads x 4


For this lab we will work only with a subset of the data (first 1000 reads): 

	head -n 4000 SRR518875.fastq > example.fastq

### Alignment using bowtie

We need the Yeast genome index files from the Bowtie website:

	cd ~
	mkdir bowtie_indexes
	cd bowtie_indexes
	curl -O ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/s_cerevisiae.ebwt.zip
	unzip s_cerevisiae.ebwt.zip

We are ready to align the reads with the default parameters of Bowtie:
	
	module load bowtie
	cd ~/statgenomics/lab_seqalign
	time bowtie --sam ~/bowtie_indexes/s_cerevisiae example.fastq example.sam

Let's take a look at the resulting SAM file containing alignment results. (SAM format specification: [http://samtools.github.io/hts-specs/SAMv1.pdf](http://samtools.github.io/hts-specs/SAMv1.pdf))

	head -n 25 example.sam

Note: to display bowtie options, you can just enter 'bowtie' at the command line to see a menu of options and the required syntax:

	bowtie

### More information on Bowtie 1 and Bowtie 2

We used Bowtie 1 in this example because the reads are only 36bp. Most reads nowadays are longer so the software authors generally recommend Bowtie 2: [http://bowtie-bio.sourceforge.net/bowtie2/faq.shtml](http://bowtie-bio.sourceforge.net/bowtie2/faq.shtml)

Bowtie 1 indexes for the human genome can be found at: [ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/](ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/). Look for hgXX files. For the most recent version of the human genome look for hg19_c.ebwt.1.zip and hg19_c.ebwt.2.zip for the paired end read indexes and hg19_c.ebwt.zip for the unpaired read index.

Bowtie 2 indexes can be found at: [ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/). Bowtie 2 indexes consist of 6 files that can be downloaded all together in one large zip archive from [ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip) or in 3 smaller parts: hg19.1.zip, hg19.2.zip, hg19.3.zip.

### Introduction to Samtools

Samtools is a nice suite of tools that allows for exploration and analysis of sequencing data.

Samtools documentation: [http://www.htslib.org/doc/samtools.html](http://www.htslib.org/doc/samtools.html)

SAM files are typically very large - particularly for large sequencing experiments in humans. It's best practice to convert SAM files to BAM files to save space:

	module load samtools
	samtools view -b example.sam -o example.bam

BAM files are binary versions of SAM files that are much smaller and can be efficiently manipulated to look at specific genomic locations. We can view alignment information from BAM files with `samtools view`:

	samtools view example.bam | head

We can compute coverage over a given region with `samtools depth`, but first we must `sort` and `index` our BAM file:

	samtools sort -o example_sorted.bam -O bam -T temp example.bam
	samtools index -b example_sorted.bam

Let's look at the coverage in a region indicated by the first read of the SAM file:

	samtools depth -r Scchr04:114400-115400 example_sorted.bam

### Running batch jobs on the cluster

In this directory, I've included a shell script that runs bowtie as we did in this lab. You can submit this as a batch job with:

	qsub run_bowtie.sh

You can check on the status of the job with `qstat`.

This is useful for much longer running jobs. For more information about running jobs on our cluster, visit [http://www.jhpce.jhu.edu/knowledge-base/how-to/](http://www.jhpce.jhu.edu/knowledge-base/how-to/). In particular, take a look at the memory requirements section: [http://www.jhpce.jhu.edu/knowledge-base/how-to/#MemSpec](http://www.jhpce.jhu.edu/knowledge-base/how-to/#MemSpec). The default limits for `mf`, `h_vmem`, and `h_fsize` are all bigger than what I specified in the shell script, but these will definitely need to be set larger for real experiments.