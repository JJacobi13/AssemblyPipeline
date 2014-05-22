"""
@author: Jetse
@version: 0.1
@attention: TODO: Remove all # from the FileControl methods to include the control (biopython is needed for testing)
@attention: Not tested yet...
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
            "minCounts": 10
        }
    }
}
"""
###############
##  Imports  ##
###############
# from qualityControl import FileControl

#############
##  Quake  ##
#############
##  Paired end
rule quakePaired:
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
        shell("correct -f {input.fileNames} -k {kmer} -c {minCounts} -m {input.countsFile} -p {threads}".format(input=input,
                                                                                                                kmer=CONFIG["kmerSize"],
                                                                                                                minCounts=CONFIG["options"]["quake"]["minCounts"],
                                                                                                                threads=threads))
        os.rename(output.forward.replace(".quake_1","_1.cor"), output.forward)
        os.rename(output.reversed.replace(".quake_2","_2.cor"), output.reversed)
#         FileControl.fastqControl(output.forward, output.reversed)


rule kmerCountsPaired:
    input:
        forward = "{samples}_1.fastq",
        reversed = "{samples}_2.fastq"
    output: "{samples}.counts"
    shell: "cat {input.forward} {input.reversed} | count-qmers -k " + CONFIG.getGlobalOption("kmer") + " > {output[0]}"

#Put paired end data file names into a file with a whitespace inbetween.
rule fastqNamesFile:
    input:
        forward = "{samples}_1.fastq",
        reversed = "{samples}_2.fastq"
    output:"{samples}.filenames.txt"
    shell: "echo \"{input.forward} {input.reversed}\" > {output[0]}"
            
##  Unpaired
rule quakeSingle:
    input:
        fastq = "{sample}.fastq",
        fileNames = "{sample}.filenames.txt",
        countsFile = "{sample}.counts"
    output: "quake.{sample}.fastq"
    threads: 999
    run:
        shell("correct -f {input.fileNames} -k {kmer} -c {minCounts} -m {input.countsFile} -p {threads}".format(input=input,
                                                                                                                kmer=CONFIG["kmerSize"],
                                                                                                                minCounts=CONFIG["options"]["quake"]["minCounts"],
                                                                                                                threads=threads))
#         FileControl.fastqControl(output.forward, output.reversed)
        
rule kmerCountsSingle:
    input: "{sample}.fastq"
    output: "{sample}.counts"
    shell: "cat {input[0]} | count-qmers -k " + CONFIG.getGlobalOption("kmer") + " > {output[0]}"       
        
#Write the filename to a file        
rule fastqNameFile:
    input: "{sample}.fastq"
    output: "{sample}.filenames.txt"
    shell: "echo \"{input[0]}\" > {output[0]}"
            
            