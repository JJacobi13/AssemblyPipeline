"""
@author: Jetse
@version: 0.2
 
Trimming of fastq files. For usage, include this in your workflow.

Required programs:
* trimmomatic
* ea-utils

Expects a global variable CONFIG (e.g. parsed from json) of at least the following structure:
{
	"illuminaAdapters": "/path/to/adapters",
	"options":{
		"fastqMcf":{
			"dupLen": "60",
			"optionalOptions": ""
		},
		
		"trimmomatic":{
			"jar": "/home/jaco001/programs/trimmomatic-0.32.jar",
			"seedMisMatches": 2,
			"palindromeClipTreshold": 30,
			"simpleClipThreshhold": 10,
			"LeadMinTrimQual": 3,
			"TrailMinTrimQual": 3,
			"windowSize": 4,
			"avgMinQual": 15,
			"minReadLen": 36,
			"phred" : "-phred33"
		}
	}
}
Expects the input files in the working directory. For paired end data the files have to end with _1.fastq (forward reads) and 
_2.fastq (reversed reads).

"""
###############
##  Imports  ##
###############
from qualityControl import FileControl
import collections

#################
##  Fastq-mcf  ##
#################
ruleorder: fastqMcfPaired > fastqMcfSingle

rule fastqMcfPaired:
	input: 
		forward = "{sample}_1.fastq",
		reversed = "{sample}_2.fastq"
	output: 
		forward = "mcf.{sample}_1.fastq",
		reversed = "mcf.{sample}_2.fastq"
	run: 
		shell("fastq-mcf {optionalParams} -D  {dups} -o {output.forward} -o {output.reversed} {adapterFile} {input.forward} {input.reversed}".format(optionalParams=CONFIG["options"]["fastqMcf"]["optionalOptions"],
																																					dups=CONFIG["options"]["fastqMcf"]["dupLen"],
																																					adapterFile=CONFIG["illuminaAdapters"],
																																					input=input,
																																					output=output))
		FileControl.fastqControl(output.forward, output.reversed)
		
rule fastqMcfSingle:
	input:"{sample}.fastq",
	output: "mcf.{sample}.fastq"
	run: 
		shell("fastq-mcf {optionalParams} -D  {dups} -o {output[0]} {adapterFile} {input[0]}".format(optionalParams=CONFIG["options"]["fastqMcf"]["optionalOptions"],
																										dups=CONFIG["options"]["fastqMcf"]["dupLen"],
																										adapterFile=CONFIG["illuminaAdapters"],
																										input=input,
																										output=output))
		FileControl.fastqControl(output[0])

###################
##  Trimmomatic  ##
###################
ruleorder: trimmomaticPaired > trimmomaticSingle

rule trimmomaticPaired:
	input:
		forward = "{sample}_1.fastq",
		reversed = "{sample}_2.fastq"
	output:
		forward = "trim.{sample}_1.fastq",
		reversed = "trim.{sample}_2.fastq",
		forwardUnpaired = "trim.{sample}_unpaired_1.fastq",
		reversedUnpaired = "trim.{sample}_unpaired_2.fastq"
	threads: 999
	run: 
		TrimOpts = collections.namedtuple("TrimOpts", CONFIG["options"]["trimmomatic"].keys())
		trimOpts = TrimOpts(**CONFIG["options"]["trimmomatic"])
		shell("java -jar {trimOpts.jar} PE {trimOpts.phred} -threads {threads} "
			"{input.forward} {input.reversed} {output.forward} {output.forwardUnpaired} {output.reversed} {output.reversedUnpaired}  "
			"ILLUMINACLIP:{adapterFile}:{trimOpts.seedMisMatches}:{trimOpts.palindromeClipTreshold}:{trimOpts.simpleClipThreshhold} "
			"LEADING:{trimOpts.LeadMinTrimQual} "
			"TRAILING:{trimOpts.TrailMinTrimQual} "
			"SLIDINGWINDOW:{trimOpts.windowSize}:{trimOpts.avgMinQual} "
			"MINLEN:{trimOpts.minReadLen}"
			"".format(
					input=input,
					output=output,
					adapterFile=CONFIG["illuminaAdapters"],
					threads=threads,
					trimOpts=trimOpts))
		FileControl.fastqControl(output.forward, output.reversed)
		
rule trimmomaticSingle:
	input: "{sample}.fastq"
	output: "trim.{sample}.fastq"
	threads: 999
	run:
		TrimOpts = collections.namedtuple("TrimOpts", CONFIG["options"]["trimmomatic"].keys())
		trimOpts = TrimOpts(**CONFIG["options"]["trimmomatic"])
		shell("java -jar {trimOpts.jar} SE {trimOpts.phred} -threads {threads} "
			"{input[0]} {output[0]} "
			"ILLUMINACLIP:{adapterFile}:{trimOpts.seedMisMatches}:{trimOpts.palindromeClipTreshold}:{trimOpts.simpleClipThreshhold} "
			"LEADING:{trimOpts.LeadMinTrimQual} "
			"TRAILING:{trimOpts.TrailMinTrimQual} "
			"SLIDINGWINDOW:{trimOpts.windowSize}:{trimOpts.avgMinQual} "
			"MINLEN:{trimOpts.minReadLen}".format(
					input=input,
					output=output,
					adapterFile=CONFIG["illuminaAdapters"],
					threads=threads,
					trimOpts=trimOpts))
		FileControl.fastqControl(output[0])
	
# rule cutadapt:
# 	input:
# 		forward = "{samples}_1.fastq",
# 		reversed = "{samples}_2.fastq"
# 	output:
# 		forward = "{samples}.cut_1.fastq",
# 		reversed = "{samples}.cut_2.fastq"
# 	shell: "touch {output.forward} {output.reversed}"


	



		



