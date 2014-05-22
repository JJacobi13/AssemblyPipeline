## Written by Jetse Jacobi
## Usage: snakemake --snakefile ~/programmingProjects/assembly/snakemake/main/assembly.py [target [target ...]]
## 
## Targets:
## wgs - wgs assembly
## allpaths - allpaths assembly
## pbJelly - reads scaffolded with PBJelly (with pacbio reads)
## scarpa - reads scaffolded with scarpa (with paired end reads)
##
## Example: snakemake --snakefile ~/programmingProjects/assembly/snakemake/main/assembly.py contigs.allpaths.fasta
## This example executes the following programs: prepareAllpathsInput.pl -> allpaths assembler
##
## Example 2: snakemake --snakefile ~/programmingProjects/assembly/snakemake/main/assembly.py wgs.pbJelly.fasta
## This example executes the following programs: fastqToCA -> wgs assembler -> PBJelly
##
################################
## Configuration requirements ##
## The following global options are required:
## outputDirector, maxThreads, expCoverage, maxMem, kmer and overwrite
##
## Each library has to be determined with [library]
## The first attribute of each library has to be libName = {library name}
## Other required attributes for each single end libraries: forward, format, sequencingPlatform, readlen
## Required attribute for paired end libraries: forward, reversed, format, type, insertSize, stdev, readlen, sequencingplatform 
##
######################################################################################################################################################
##  Imports  ##
###############




################
##  Allpaths  ##
################
rule prepareAllpaths:
    input: dynamic("{sample}.fastq")
    output: "somepath/something.fastb"
    shell: "touch {output[0]}"  

rule allpaths:
    input: rules.prepareAllpaths.output
    output: "assembly/assembly/assembly/ASSEMBLIES/assembly/contigs.fasta"
    shell: "mkdir -p assembly/assembly/assembly/ASSEMBLIES/assembly/; touch {output[0]}"
 
rule allpathsCleanup:
    input: rules.allpaths.output
    output: 
        contigs = "allpaths.contigs.fasta",
        scaffolds = "allpaths.scaffolds.fasta"
    shell: "rm -r assembly; touch {output.contigs}; touch {output.scaffolds}"

################
##  PB jelly  ##
################
rule pbJelly:
    input:
        assembly = "{assembler}.contigs.fasta",
        pacbio = "pacbio.fastq"
    output: "{assembler}.pbJelly.scaffolds.fasta"
    shell: "touch {output[0]}"

##############
##  Scarpa  ##
##############
#Strange output format because all output files have to contain all wildcards...
rule scarpaProcess:
    input:
        assembly = "{assembler}.contigs.fasta",
        forward="{sample}_1.fastq",
        reversed="{sample}_2.fastq"
    output: 
        fasta="{assembler}.contigs.{sample}.scarpa.fa",
        info="{assembler}.contigs.{sample}.scarpa.info",
        forward="{assembler}.{sample}_1.fastq.scarpa.fq",
        reversed="{assembler}.{sample}_2.fastq.scarpa.fq"
    shell: "touch {output.fasta}; touch {output.info}; touch {output.forward}; touch {output.reversed}"

rule scarpaParse:
    input:
        assembly = "{assembler}.contigs.{sample}.scarpa.fa",
        forward="{assembler}.{sample}_1.fastq.scarpa.fq",
        reversed="{assembler}.{sample}_2.fastq.scarpa.fq"
    output: "{assembler}.{sample}.map"
    shell: "touch {output[0]}"
    
rule scarpa:
    input:
        assembly="{assembler}.contigs.{sample}.scarpa.fa",
        map="{assembler}.{sample}.map",
        info="{assembler}.contigs.{sample}.scarpa.info"
    output: "{assembler}.{sample}.scarpa.scaffolds.fasta"
    shell: "touch {output[0]}"
     


    