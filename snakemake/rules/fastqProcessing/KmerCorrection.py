"""
@author: Jetse
@version: 0.2
@attention: Not tested yet, full dataset needed...
@bug: output of quake single probably has to be changed, will figure out in tests

"Quake is a package to correct substitution sequencing errors in experiments with deep coverage (e.g. >15X), 
specifically intended for Illumina sequencing reads. Quake adopts the k-mer error correction framework, 
first introduced by the EULER genome assembly package. Unlike EULER and similar progams, 
Quake utilizes a robust mixture model of erroneous and genuine k-mer distributions to determine where errors are 
located. Then Quake uses read quality values and learns the nucleotide to nucleotide error rates to determine what 
types of errors are most likely. This leads to more corrections and greater accuracy, especially with respect to 
avoiding mis-corrections, which create false sequence unsimilar to anything in the original genome sequence from 
which the read was taken."

Expects a global variable CONFIG (e.g. parsed from json) of at least the following structure:
{
    "kmerSize": 17
    "options":{
        "quake":{
            "path": "/home/VLPB/programs/Quake/bin/quake.py",
            "optionalOpts": ""
        }
    }
}
"""
###############
##  Imports  ##
###############
from qualityControl import FileControl

#############
##  Quake  ##
#############
#If data is paired end, use the paired end rules...
ruleorder: quakePaired > quakeSingle

##  Paired end
rule quakePaired:
    input: "{sample}.filenames.txt",
    output:
        forward = "quake.{sample}_1.fastq",
        reversed = "quake.{sample}_2.fastq" 
    threads: 999
    run: 
        shell("python2 {quakePath} {optional} -f {readsFile} -k {kmerSize} -p {threads}"
              "".format(quakePath=CONFIG["options"]["quake"]["path"],
                        optional=CONFIG["options"]["quake"]["optionalOpts"],
                        readsFile=input[0],
                        kmerSize=CONFIG["kmerSize"],
                        threads=threads))
        os.rename(wildcards.sample + "_1.cor.fastq", output.forward)
        os.rename(wildcards.sample + "_2.cor.fastq", output.reversed)
        FileControl.fastqControl(output.forward, output.reversed)

#Put paired end data file names into a file with a whitespace inbetween.
rule fastqNamesFile:
    input:
        forward = "{samples}_1.fastq",
        reversed = "{samples}_2.fastq"
    output:"{samples}.filenames.txt"
    shell: "echo \"{input.forward} {input.reversed}\" > {output[0]}"
            
##  Unpaired
rule quakeSingle:
    input: fastq = "{sample}.fastq"
    output: "quake.{sample}.fastq"
    threads: 999
    run:
        shell("python2 {quakePath} {optional} -r {readsFile} -k {kmerSize} -p {threads}"
              "".format(quakePath=CONFIG["options"]["quake"]["path"],
                        optional=CONFIG["options"]["quake"]["optionalOpts"],
                        readsFile=input[0],
                        kmerSize=CONFIG["kmerSize"],
                        threads=threads))
        os.rename(wildcards.sample + ".cor.fastq", output[0])
        FileControl.fastqControl(output[0])