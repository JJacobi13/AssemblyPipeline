"""
@author: Jetse
@version: 0.1
@attention: TODO: Remove all # from the FileControl methods to include the control (biopython is needed for testing)
@attention: Not tested yet...

This module changes all accepted inputfiles -- determined by their suffixes -- to fastq format.
All accepted input formats are:
+ *.sff -- 454 raw input files
+ *.fq -- fastq files, only a symlink to the original file is created

"""

#################################
##  fq suffix to fastq suffix  ##
#################################
rule fqToFastq:
    input: "{prefix}.fq"
    output: "{prefix}.fastq"
    shell: "ln -s {input[0]} {output[0]}"
    
####################
##  SFF to fastq  ##
####################
rule sffToFastq:
    input: inFile = "{samples}.sff"
    output: outFile = "{samples}.fastq"
    shell: "sff2fastq -o {output.outFile} {input.inFile}"
    