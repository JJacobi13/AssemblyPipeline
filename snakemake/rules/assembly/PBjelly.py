"""
@author: Jetse
@version: 0.1
@attention: Not tested yet...

"PBJelly is a highly automated pipeline that aligns long sequencing reads (such as PacBio RS reads or long 454 reads in fasta format) 
to high-confidence draft assembles. PBJelly fills or reduces as many captured gaps as possible to produce upgraded draft genomes."

 
"""
################
##  PB jelly  ##
################
rule pbJelly:
    input:
        assembly = "{assembler}.contigs.fasta",
        pacbio = "pacbio.fastq"
    output: "{assembler}.pbJelly.scaffolds.fasta"
    shell: "touch {output[0]}"
