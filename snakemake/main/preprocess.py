## Written by Jetse Jacobi
## Usage: snakemake --snakefile ~/programmingProjects/assembly/snakemake/main/preprocess.py [target [target ...]]
## 
## Targets:
## + mcf - Trimmed with fastq-mcf
## + {organism}filtered - {organism} reads are removed with bowtie2, the reference fasta has to be in REF_GENOMES
## + quake - Quake corrected reads (only works for reads shorter than 125 bp)
## + merged - Fragment reads merged with SeqPrep
## + trim - Trimmomatic trimmed reads
## - cut - Cutadapt adapter removed
##
## Example: snakemake --snakefile ~/programmingProjects/assembly/snakemake/main/preprocess.py frag180.mcf.quake.merged.fastq 
## This example executes the following programs: fastq-mcf -> quake correction -> seqprep
##
## Example 2: snakemake --snakefile ~/programmingProjects/assembly/snakemake/main/preprocess.py frag180.quake.trim_1.fastq frag180.quake.trim_2.fastq 
## This example executes the following programs: quake correction -> trimmomatic on the paired end reads
##
## Configuration requirements:
## File in working directory, name: config.ini
##
## File requirements:
## Before first dot has to be the library name defined in the configuration file
## Paired end reads have to end with _1.fastq and _2.fastq
## 
#################################################################################################################################
##  Imports  ##
###############
from configuration import Configuration
from pipelineUtils import FastqUtils
import os
from qualityControl import FileControl
################
##  Settings  ##
################
CONFIG = Configuration.instance
CONFIG.parseIni("config.ini")

MIN_READ_LEN = "36"
DUPLICATE_REQ_OVERLAP = "60"

#trimming
MIN_TRIM_QUALITY="3"

#clipping
ADAPTER_MISMATCHES = "2"

#Quake specific
QUAKE_MINIMUM_COUNTS = "20"

######################
##  Required files  ##
######################
ADAPTERFILE = "/home/jaco001/data/illumina-adapter-trimming-file.fna"
REF_GENOMES = {"PhiX":"/home/jaco001/programmingProjects/assembly/src/preprocessing/phixDb/PhiX.fasta"} 
TRIMMOMATIC_JAR = "/home/jaco001/programs/trimmomatic-0.32.jar"
TRIMMOMATIC_PHRED_ENCODING = "-phred33"

###################################################################
##  Adapter clipping / bad quality trimming / duplicate removal  ##
###################################################################
rule fastqMcf:
	input: 
		forward = "{samples}_1.fastq",
		reversed = "{samples}_2.fastq"
	output: 
		forward = "{samples}.mcf_1.fastq",
		reversed = "{samples}.mcf_2.fastq"
	run: 
		shell("fastq-mcf -D  {DUPLICATE_REQ_OVERLAP} -o {output.forward} -o {output.reversed} {ADAPTERFILE} {input.forward} {input.reversed}")
		FileControl.fastqControlPaired(output.forward, output.reversed)

rule trimmomatic:
	input:
		forward = "{samples}_1.fastq",
		reversed = "{samples}_2.fastq"
	output:
		forward = "{samples}.trim_1.fastq",
		reversed = "{samples}.trim_2.fastq",
		forwardUnpaired = "{samples}.trim_unpaired_1.fastq",
		reversedUnpaired = "{samples}.trim_unpaired_2.fastq",
	threads: 5
	run: 
		shell("java -jar {TRIMMOMATIC_JAR} PE {TRIMMOMATIC_PHRED_ENCODING} -threads {threads} {input.forward} {input.reversed} {output.forward} {output.forwardUnpaired} {output.reversed} {output.reversedUnpaired}  ILLUMINACLIP:{ADAPTERFILE}:{ADAPTER_MISMATCHES}:30:10 LEADING:{MIN_TRIM_QUALITY} TRAILING:{MIN_TRIM_QUALITY} SLIDINGWINDOW:4:15 MINLEN:{MIN_READ_LEN}")
		FileControl.fastqControlPaired(output.forward, output.reversed)
	
# rule cutadapt:
# 	input:
# 		forward = "{samples}_1.fastq",
# 		reversed = "{samples}_2.fastq"
# 	output:
# 		forward = "{samples}.cut_1.fastq",
# 		reversed = "{samples}.cut_2.fastq"
# 	shell: "touch {output.forward} {output.reversed}"

########################
##  File conversions  ##
########################
rule sffToFastq:
	input: inFile = "{samples}.sff"
	output: outFile = "{samples}.fastq"
	shell: "sff2fastq -o {output.outFile} {input.inFile}"
	
###############################
##  Contamination filtering  ##
###############################	
rule filterPhiX:
	input: 
		forward = "{samples}_1.fastq",
		reversed = "{samples}_2.fastq",
		reference = lambda wildcards: REF_GENOMES[wildcards.org],
		index = lambda wildcards: REF_GENOMES[wildcards.org] + ".1.bt2"
	output: 
		forward = "{samples}.{org}filtered_1.fastq",
		reversed = "{samples}.{org}filtered_2.fastq"
	threads: 999
	run: 
		insertSize = int(CONFIG.getLibInfo(input.forward.split(".")[0][:-2]).getOption("insertSize"))
		minInsert = str(insertSize*0.5)
		maxInsert = str(insertSize*2)
		shell("bowtie2 -p {threads} -x {input.reference} --un-conc {output.forward}_tmp_unmapppedPhix.fastq -I {minInsert} -X {maxInsert} -1 {input.forward} -2 {input.reversed} -S /dev/null")
		os.rename(output.forward + "_tmp_unmapppedPhix.1.fastq", output.forward)
		os.rename(output.forward + "_tmp_unmapppedPhix.2.fastq", output.reversed)
		FileControl.fastqControlPaired(output.forward, output.reversed)
		
	
rule bowtieIndex:
	input: inFile = "{ref}.fasta"
	output: outFile = "{ref}.fasta.1.bt2"
	shell: "bowtie2-build {input.inFile} {input.inFile}"

#################################
##  Merging overlapping reads  ##
#################################
rule seqprep:
	input:
		forward = "{samples}_1.fastq",
		reversed = "{samples}_2.fastq"
	output: 
		merged = "{samples}.merged.fastq",
		forwardSingle = "{samples}.seqSingle_1.fastq",
		reversedSingle = "{samples}.seqSingle_2.fastq"
	run: 
		phredEncoding = "" if FastqUtils.determineQuality(input.forward) == 32 else "-6 "
		shell("SeqPrep " + phredEncoding + "-f " + input.forward + " -r " + input.reversed + " -1 "+output.forwardSingle+".gz -2 " + output.reversedSingle + ".gz -s " + output.merged + ".gz")
		#TODO: Find better way to do unzipping with the rule, and still execute fastq control on all output files...
		shell("gunzip " + output.merged + ".gz")
		shell("gunzip " + output.forwardSingle + ".gz")
		shell("gunzip " + output.reversedSingle + ".gz")
		FileControl.fastqControl(output.merged)
		FileControl.fastqControl(forwardSingle.merged)
		FileControl.fastqControl(reversedSingle.merged)
		
#######################
##  kmer correction  ##
#######################
rule quake:
	input:
		forward = "{samples}_1.fastq",
		reversed = "{samples}_2.fastq",
		fileNames = "{samples}.filenames.txt",
		countsFile = "{samples}.counts"
	output:
		forward = "{samples}.quake_1.fastq",
		reversed = "{samples}.quake_2.fastq" 
	threads: 999
	run: 
		shell("correct -f {input.fileNames} -k " + CONFIG.getGlobalOption("kmer") + " -c {QUAKE_MINIMUM_COUNTS} -m {input.countsFile} -p {threads}")
		os.rename(output.forward.replace(".quake_1","_1.cor"), output.forward)
		os.rename(output.reversed.replace(".quake_2","_2.cor"), output.reversed)
		FileControl.fastqControlPaired(output.forward, output.reversed)

rule kmerCounts:
	input:
		forward = "{samples}_1.fastq",
		reversed = "{samples}_2.fastq"
	output: "{samples}.counts"
	shell: "cat {input.forward} {input.reversed} | count-qmers -k " + CONFIG.getGlobalOption("kmer") + " > {output[0]}"
	
###################
##  Other rules  ##
###################	
#Put paired end data file names into a file with a whitespace inbetween.
rule fastqNamesFile:
	input:
		forward = "{samples}_1.fastq",
		reversed = "{samples}_2.fastq"
	output:"{samples}.filenames.txt"
	run:
		with open(output[0], "w") as writer:
			writer.write(input.forward + " ")
			writer.write(input.reversed)
#First with .fastq.gz, because .gz gives maximum recursion depth exceeded... 
#TODO: find better way, hate workarounds...
rule unzip:
	input: "{file}.fastq.gz"
	output: "{file}.fastq"
	shell: "gunzip {input[0]}"

