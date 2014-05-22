import json, os
 
with open("config.json") as f:
    CONFIG = json.load(f)

## Create a symlink so the data is perfectly fit for all rules. All reads contain the name of their library. 
## All forward reads are {library name}_1.fastq and all reversed reads are {library name}_2.fastq
for sample in CONFIG["libraries"]:
    for i in range(len(CONFIG["libraries"][sample]["reads"])):
        try:
            os.symlink(CONFIG["libraries"][sample]["reads"][i], sample +"_" + str(i+1) + ".fastq")
        except FileExistsError:
            pass

#TODO: relative paths, __file__ is path to snakemake...       
include: "/home/jaco001/programmingProjects/assembly/snakemake/rules/fastqProcessing/Trimming.py"
 
include: "/home/jaco001/programmingProjects/assembly/snakemake/rules/fastqProcessing/ContaminationFiltering.py"

