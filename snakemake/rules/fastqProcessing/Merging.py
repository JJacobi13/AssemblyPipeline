'''
Created on May 16, 2014

@author: Jetse
'''
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