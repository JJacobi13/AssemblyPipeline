The fastq processing package
============================
The fastq processing package contains rules for tools which are able to improve the quality of fastq reads. 
Also rules for visualising the quality of fastq reads can be found in this package.

Prefixes
----------
Each rule is determined whether it has to be executed by prefixes in the filename. 
All prefixes are shown in the following list:  

  * Trimming module
  	* mcf - Trim fastq reads with fastq-mcf
  	* trim - Trim fastq reads with Trimmomatic
  * ContaminationFiltering module (untested)
  	* {organism}Filtered - {organism} reads are removed with bowtie2, more info below in section Organism filtering
  * KmerCorrection module
  	* quake - Correct fastq reads with Quake
  * Merging module
  	* merged - Fragment fastq reads merged with SeqPrep
  	
Suffixes
--------
To make all rules work, all data has to be converted to fastq format. All supported suffixes are:

  * .sff (454 reads) are converted with sff2fastq
  * .fq - a symlink is created to this file
   
Fastq-mcf
---------
> Fastq-mcf detects levels of adapter presence, computes likelihoods and locations (start, end) of the adapters. 
Removes the adapter sequences from the fastq file(s). fastq-mcf attempts to:

> * Detect & remove sequencing adapters and primers
* Detect limited skewing at the ends of reads and clip
* Detect poor quality at the ends of reads and clip
* Detect Ns, and remove from ends
* Remove reads with CASAVA 'Y' flag (purity filtering)
* Discard sequences that are too short after all of the above
* Keep multiple mate-reads in sync while doing all of the above

Organism filtering
------------------
The organism filtering can be executed for any organism with a known reference genome. 
The {organism} is used as a wildcard, and has to be present as a key in the JSON file below the 
"contaminationRefGenomes". A small example is shown below:

	"contaminationRefGenomes" : {
		"PhiX" : "/home/jaco001/programmingProjects/assembly/src/preprocessing/phixDb/PhiX.fasta"
		"eColi" : "/home/jaco001/programmingProjects/assembly/src/preprocessing/eColiDb/eColi.fasta"
	},
	
In this JSON there are two reference genomes. When the fastq reads which map on the PhiX reference genome have to be
filtered, the prefix PhiXFiltered has to be used. The command eColiFiltered will filter all reads out that map on the
eColi reference genome.

Quake
-----
> Quake is a package to correct substitution sequencing errors in experiments with deep coverage (e.g. >15X), 
specifically intended for Illumina sequencing reads. Quake adopts the k-mer error correction framework, 
first introduced by the EULER genome assembly package. Unlike EULER and similar progams, 
Quake utilizes a robust mixture model of erroneous and genuine k-mer distributions to determine where errors are 
located. Then Quake uses read quality values and learns the nucleotide to nucleotide error rates to determine what 
types of errors are most likely. This leads to more corrections and greater accuracy, especially with respect to 
avoiding mis-corrections, which create false sequence unsimilar to anything in the original genome sequence from 
which the read was taken.

SeqPrep
-------
> SeqPrep is a program to merge paired end Illumina reads that are overlapping into a single longer 
read. It may also just be used for its adapter trimming feature without doing any paired end overlap. 
When an adapter sequence is present, that means that the two reads must overlap (in most cases) 
so they are forcefully merged. When reads do not have adapter sequence they must be treated with 
care when doing the merging, so a much more sensitive approach is taken. The default parameters 
were chosen with sensitivity in mind, so that they could be ran on libraries where very few reads are 
expected to overlap. It is always safest though to save the overlapping procedure for libraries where 
you have some prior knowledge that a significant portion of the reads will have some overlap.

Trimmomatic
-----------
> Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.
The selection of trimming steps and their associated parameters are supplied on the command line.
The current trimming steps are:

> * ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
* SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
* LEADING: Cut bases off the start of a read, if below a threshold quality
* TRAILING: Cut bases off the end of a read, if below a threshold quality
* CROP: Cut the read to a specified length
* HEADCROP: Cut the specified number of bases from the start of the read
* MINLEN: Drop the read if it is below a specified length
* TOPHRED33: Convert quality scores to Phred-33
* TOPHRED64: Convert quality scores to Phred-64
 
> It works with FASTQ (using phred + 33 or phred + 64 quality scores, depending on the Illumina pipeline used), 
either uncompressed or gzipp'ed FASTQ. Use of gzip format is determined based on the .gz extension.

> For single-ended data, one input and one output file are specified, plus the processing steps. 
For paired-end data, two input files are specified, and 4 output files, 2 for the 'paired' output where both 
reads survived the processing, and 2 for corresponding 'unpaired' output where a read survived, but the partner 
read did not.
