"""
@author: Jetse
@version: 0.1
@attention: Not tested yet...

ALLPATHS-LG is a short read assembler and it works on both small and large (mammalian size) genomes. 
To use it, you should first generate ~100 base Illumina reads from two libraries: 
one with ~180 bp insert size (fragment data), and one with a ~3000 bp insert size (jumping library), both at about 45x coverage. 
Sequence from longer fragments will enable longer-range continuity.



"""

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
