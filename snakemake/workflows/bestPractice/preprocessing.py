import json, os, sys
 
with open("config.json") as f:
    CONFIG = json.load(f)

## Create a symlink so the data is perfectly fit for all rules. All reads contain the name of their library. 
## All forward reads are {library name}_1.fastq and all reversed reads are {library name}_2.fastq
for sample in CONFIG["libraries"]:
    for i in range(len(CONFIG["libraries"][sample]["reads"])):
        try:
            if os.path.exists("preprocessing") == False:
                os.mkdir("preprocessing")
            os.symlink(CONFIG["libraries"][sample]["reads"][i], "preprocessing/"+ sample +"_" + str(i+1) + ".fastq")
        except FileExistsError:
            pass

#TODO: relative paths, __file__ is path to snakemake...       
include: os.path.dirname(sys.argv[2])+"/../../rules/fastqProcessing/Trimming.py"
 
include: os.path.dirname(sys.argv[2])+"/../../rules/fastqProcessing/ContaminationFiltering.py"

include: os.path.dirname(sys.argv[2])+"/../../rules/fastqProcessing/Merging.py"

include: os.path.dirname(sys.argv[2])+"/../../rules/fastqProcessing/KmerCorrection.py"

#Assembly, preprocessing is wrong name :P
include: os.path.dirname(sys.argv[2])+"/../../rules/assembly/Wgs.py"

rule all:
    input: "assembly/wgs.contigs.fasta"