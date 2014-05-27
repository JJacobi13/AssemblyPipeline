"""
@author: Jetse
@version: 0.2
 
Merging of fragment data. For usage, include this in your workflow.

Required programs:
* SeqPrep

Expects a global variable CONFIG (e.g. parsed from json) of at least the following structure:
{
    "options":{
        "SeqPrep":{
            "phredEncoding": "-6", #"-6" for phred 64, "" for phred 33
            "optionalOptions": ""
        }
    }
}
Expects the input files in the working directory. The input files have to end with "_1" 
(forward reads) and "_2" (reversed reads) 

"""
#################################
##  Merging overlapping reads  ##
#################################
rule seqprep:
    input:
        forward = "preprocessing/{samples}_1.fastq",
        reversed = "preprocessing/{samples}_2.fastq"
    output: 
        merged = "preprocessing/merged.{samples}.fastq",
        forwardSingle = "preprocessing/seqSingle.{samples}_1.fastq",
        reversedSingle = "preprocessing/seqSingle.{samples}_2.fastq"
    run: 
        shell("SeqPrep {phredEncoding} {optional}"
              "-f {input.forward} -r {input.reversed} "
              "-1 {output.forwardSingle}.gz -2 {output.reversedSingle}.gz "
              "-s {output.merged}.gz".format(phredEncoding=CONFIG["options"]["SeqPrep"]["phredEncoding"],
                                             optional=CONFIG["options"]["SeqPrep"]["optionalOptions"],
                                             input=input,
                                             output=output))
        #TODO: Find better way to do unzipping with the rule, and still execute fastq control on all output files...
        shell("gunzip " + output.merged + ".gz")
        shell("gunzip " + output.forwardSingle + ".gz")
        shell("gunzip " + output.reversedSingle + ".gz")
        FileControl.fastqControl(output.merged)
        FileControl.fastqControl(output.forwardSingle)
        FileControl.fastqControl(output.reversedSingle)