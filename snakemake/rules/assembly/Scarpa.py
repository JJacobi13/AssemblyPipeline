"""
@author: Jetse
@version: 0.1
@attention: Not tested yet...

"Scarpa is a stand-alone scaffolding tool for NGS data with paired end reads. It can be used together with virtually any genome assembler and any NGS read mapper 
that supports SAM format. Other features include support for multiple libraries and an option to estimate insert size distributions from data. 
Scarpa is available free of charge for academic and commercial use under the GNU General Public License (GPL)."

{
    "libraries"{
        "libName":{
            "insertSize":500
        }
    }
}
"""
###############
##  Imports  ##
###############
import os

##############
##  Scarpa  ##
##############
#Strange output format because all output files have to contain all wildcards...
rule scarpaProcess:
    input:
        assembly = "{assembler}.contigs.fasta",
        reads="scarpaMerged.{sample}.fastq",
    output: 
        fasta="{assembler}.contigs.{sample}.scarpa.fa",
        info="{assembler}.contigs.{sample}.scarpa.info",
        reads="{assembler}.{sample}.fastq.scarpa.fq",
    run: 
    shell("scarpa_process -c {input.assembly} -f {input.reads} -i {insertSize}".format(input=input, output=output, insertSize=CONFIG["libraries"][wildcards.sample]["insertSize"]))
    os.rename(input.reads + ".scarpa.fq", output.reads)
    
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
     
rule mergeForwardReversedFastq:
    input:
        forward="{sample}_1.fastq",
        reversed="{sample}_2.fastq"
    output: "scarpaMerged.{sample}.fastq"
    run:
        revHash = []
        with open(input.forward) as ffh, open(input.reversed) as rfh:
            with open(output[0], "w") as writer:
                for fr, rr in izip(ffh, rfh):
                    writer.write(fr)
                    revHash.append(rr)
                    if len(revHash) == 4:
                        for line in revHash:
                            writer.write(line)
                            revHash = []