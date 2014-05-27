"""
@author: Jetse
@version: 0.1
@attention: TODO: Remove all # from the FileControl methods to include the control (biopython is needed for testing)
@attention: Not tested yet, bowtie2 needs 64 bit machine, virtual machine is 32 bit...

This module can be executed for any organism with a known reference genome. 
The {organism} is used as a wildcard, and has to be present as a key in the JSON file below the 
"contaminationRefGenomes". A small example is shown below:

    "contaminationRefGenomes" : {
        "PhiX" : "/home/jaco001/programmingProjects/assembly/src/preprocessing/phixDb/PhiX.fasta"
        "eColi" : "/home/jaco001/programmingProjects/assembly/src/preprocessing/eColiDb/eColi.fasta"
    },
    
In this JSON there are two reference genomes. When the fastq reads which map on the PhiX reference genome have to be
filtered, the prefix PhiXFiltered has to be used. The command eColiFiltered will filter all reads out that map on the
eColi reference genome.

Required programs:
* bowtie2

Expects a global variable CONFIG (e.g. parsed from json) of at least the following structure for paired end data:
{
    "libraries"{
        "testLibrary":{
            "insertSize":500
        }
    }
    "contaminationRefGenomes" : {
        "PhiX" : "/home/jaco001/programmingProjects/assembly/src/preprocessing/phixDb/PhiX.fasta"
    }
}
For unpaired data, the libraries section is not needed in this module.

"""
###############
##  Imports  ##
###############
# from qualityControl import FileControl

###############################
##  Contamination filtering  ##
###############################    
rule filterPaired:
    input: 
        forward = "preprocessing/{sample}_1.fastq",
        reversed = "preprocessing/{sample}_2.fastq",
        reference = lambda wildcards: CONFIG["contaminationRefGenomes"][wildcards.org],
        index = lambda wildcards: CONFIG["contaminationRefGenomes"][wildcards.org] + ".1.bt2"
    output: 
        forward = "preprocessing/{org}Filtered.{sample}_1.fastq",
        reversed = "preprocessing/{org}Filtered.{sample}_2.fastq"
    threads: 999
    run: 
        insertSize = int(CONFIG["libraries"][wildcards.sample]["insertSize"])
        print("bowtie2 -p {threads} -x {input.reference} "
              "--un-conc {output.forward}_tmp_unmapppedPhix.fastq "
              "-I {minInsert} -X {maxInsert} "
              "-1 {input.forward} -2 {input.reversed} -S /dev/null".format(threads=threads,
                                                                           input=input,
                                                                           output=output,
                                                                           minInsert=insertSize*0.5,
                                                                           maxInsert=insertSize*2))
        os.rename(output.forward + "_tmp_unmapppedPhix.1.fastq", output.forward)
        os.rename(output.forward + "_tmp_unmapppedPhix.2.fastq", output.reversed)
#         FileControl.fastqControl(output.forward, output.reversed)
            
rule filterSingle:
    input:
        fastq = "preprocessing/{sample}.fastq",
        reference = lambda wildcards: CONFIG["contaminationRefGenomes"][wildcards.org],
        index = lambda wildcards: CONFIG["contaminationRefGenomes"][wildcards.org] + ".1.bt2"
    output: "preprocessing/{org}_filtered.{sample}.fastq"
    run:
        print("bowtie2 -p {threads} -x {input.reference} "
              "--un {output[0]} "
              "-I {minInsert} -X {maxInsert} "
              "-U {input.fastq} -S /dev/null".format(threads=threads,
                                                   input=input,
                                                   output=output))
#         FileControl.fastqControl(output[0])
        
rule bowtieIndex:
    input: inFile = "{ref}.fasta"
    output: outFile = "{ref}.fasta.1.bt2"
    shell: "bowtie2-build {input.inFile} {input.inFile}"