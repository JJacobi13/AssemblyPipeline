"""
@author: Jetse
@version: 0.1
@attention: Not tested yet...

"Celera Assembler is a de novo whole-genome shotgun (WGS) DNA sequence assembler. 
It reconstructs long sequences of genomic DNA from fragmentary data produced by whole-genome shotgun sequencing. 
Celera Assembler was developed at Celera Genomics starting in 1999. It was released to SourceForge in 2004 as the 
wgs-assembler under the GNU General Public License. The pipeline revised for 454 data was named CABOG."

Expects a global variable CONFIG (e.g. parsed from json) of at least the following structure:
{
    "maxMem" : 200000000,
    "libraries":{
        "fragData":{
            "reads" : ["/home/jaco001/beauveriaBassiana/1520/raw/1520.L180_SZABPI034205-52_1.fq","/home/jaco001/beauveriaBassiana/1520/raw/1520.L180_SZABPI034205-52_2.fq"],
            "type" : "pe",
            "insertSize" : 180,
            "insertSizeStdev" : 36,
            "readlen" : 100,
            "platform" : "illumina"
        }
    },
    "wgs":{
        "input":{
            "fragData":["../preprocessing/trim.fragData_1.fastq"]
        }
    }
}

"""
###############
##  Imports  ##
###############
import math, os

##########################
##  WGS assembler rules ##
##########################
rule wgsExecution:
    input: expand("assembly/{sample}.frg",sample=CONFIG["options"]["wgs"]["input"].keys())
    output: "assembly/wgs/9-terminator/assembly.ctg.fasta"
    run: 
        specFile = SpecFile("assembly/wgs.cnf",CONFIG["maxMem"],threads)
        print("runCA -d assembly/wgs/ -p assembly -s {specFile} {inputFiles}".format(specFile=specFile.fileName, inputFiles=" ".join(input)))
        if not os.path.exists("assembly/wgs/9-terminator"):
            os.mkdir("assembly/wgs/9-terminator")
        open(output[0], "a").close()
        
        
ruleorder: fastqToCaPaired > fastqToCaSingle

rule fastqToCaPaired:
    input: lambda wildcards: CONFIG["options"]["wgs"]["input"][wildcards.sample]
    output: "assembly/{sample}.frg"
    run:
        print("fastqToCA -insertSize {insSize} {insDev} -libraryname {sample} -technology {tech} {direction} -mates {input[0]},{input[1]} "
                  "> {output[0]}".format(
                                         insSize=CONFIG["libraries"][wildcards.sample]["insertSize"],
                                         insDev=CONFIG["libraries"][wildcards.sample]["insertSizeStdev"],
                                         sample=wildcards.sample,
                                         tech=getTech(CONFIG["libraries"][wildcards.sample]["platform"], CONFIG["libraries"][wildcards.sample]["readlen"]),
                                         direction=getOrientation(CONFIG["libraries"][wildcards.sample]["type"]),
                                         input=input,
                                         output=output))
        open(output[0], 'a').close()

rule fastqToCaSingle:
    input: "preprocessing/{sample}_1.fastq"
    output: "assembly/{sample}.frg"
    run:
        print("fastqToCA -technology {tech} -reads {input[0]} > {output[0]}"
                  "".format(
                            tech=getTech(CONFIG["libraries"][wildcards.sample]["platform"], CONFIG["libraries"][wildcards.sample]["readlen"]),
                            input=input,
                            output=output
                            ))
        open(output[0], 'a').close()
      
rule wgsCleanup:
    input: "assembly/wgs/9-terminator/assembly.ctg.fasta"
    output: "assembly/wgs.contigs.fasta"
    shell: "mv {input[0]} {output[0]}"

#################
##  Functions  ##
#################
def getTech(platform, readlen):
    """
    The method getType returns the technology of the reads in the format WGS accepts.
    """
        
    if platform == "illumina":
        if readlen < 160:
            return "illumina"  
        return "illumina-long"
    elif platform == "454":
        return "454"
    elif platform == "pacbio":
        return "pacbio-raw"
    raise("Unrecognized combination of readlen and platform")
    
def getOrientation(libraryType):
        """
        The method getOrientation returns the orientation of the reads in the format WGS accepts.
        """
        if libraryType == "mp":
            return "-outtie"
        elif libraryType == "pe":
            return "-innie"
        raise ValueError("Unrecognized library type")

###############
##  Objects  ##
###############
class SpecFile(object):
    """
    The specfile object creates a specfile for the wgs assembler based on the number of threads and threads the user wants to use
    """
    def __init__(self, fileName, maxMem, maxThr):
        self.fileName = fileName
        with open(self.fileName, "w") as specWriter:
            specWriter.write("#Default allowed error rates for Illumina\n")
            specWriter.write("utgErrorRate = 0.03\n")
            specWriter.write("utgErrorLimit = 6.5\n")
            specWriter.write("cnsErrorRate = 0.06\n")
            specWriter.write("cgwErrorRate = 0.10\n")
            specWriter.write("ovlErrorRate = 0.06\n")
            specWriter.write("merSize=22\n")
            specWriter.write("unitigger=bogart\n")
            
            specWriter.write("#Don't use grid\n")
            specWriter.write("useGrid = 0\n")
            specWriter.write("scriptOnGrid = 0\n")
            specWriter.write("frgCorrOnGrid = 0\n")
            specWriter.write("ovlCorrOnGrid = 0\n")
            
            specWriter.write("#Memory settings for the given maximum memory/threads\n")
            specWriter.write("merylMemory = " + str(int(maxMem/1000)) + "\n") #merylMemory is in mb
            specWriter.write("merylThreads = " + str(maxThr) + "\n")
            
            specWriter.write("ovlStoreMemory = " + str(int(maxMem/1000)) + "\n") #ovlStoreMemory is in mb
            
            hashBits = self.getOvlHashBits(maxMem*0.66, maxThr)
            specWriter.write("ovlThreads = " + str(self.ovlThreads) + "\n")
            specWriter.write("ovlConcurrency = " + str(self.ovlConcurrency) + "\n")
            specWriter.write("ovlHashBits = " + str(hashBits) + "\n")

            specWriter.write("ovlHashBlockLength = " + str(maxMem*0.33/self.ovlConcurrency*100) + "\n")
            specWriter.write("ovlRefBlockSize = 20000000\n")
            
            specWriter.write("ovlCorrBatchSize = 1000000\n")
            specWriter.write("ovlCorrConcurrency = "+ str(int(math.floor(min(maxThr,int(math.floor(maxMem/1000))))))+"\n")
            
            frgCorrThreads = self.ovlThreads
            specWriter.write("frgCorrThreads = " + str(frgCorrThreads) + "\n")
            frgCorrConcurrency = self.ovlConcurrency
            specWriter.write("frgCorrConcurrency = " + str(frgCorrConcurrency) + "\n")
            frgCorrBatchSize = (100/float(1300))*maxMem
            specWriter.write("frgCorrBatchSize = " + str(int(frgCorrBatchSize)) + "\n")
            
            specWriter.write("merOverlapperThreads = " + str(maxThr) + "\n")
            #edit merOverlapperSeedBatchSize when documentation available
            specWriter.write("cnsConcurrency = " + str(maxThr) + "\n")
            
            specWriter.write("mbtThreads = " + str(frgCorrThreads) + "\n")
            specWriter.write("mbtConcurrency = " + str(frgCorrConcurrency) + "\n")
            
    def getOvlHashBits(self, maxMem,maxThr):
        self.ovlThreads = int(math.floor(math.sqrt(maxThr)))
        self.ovlConcurrency = int(math.floor(maxThr/self.ovlThreads))
        ovlHashBitsDict = {54:18,108:19,216:20,432:21,864:22,1728:23,3456:24,6912:25,13824:26,27648:27,55296:28,110592:29,221184:30}
        maxBitsSize = maxMem/self.ovlConcurrency/1000
        for size in sorted(ovlHashBitsDict, reverse=True):
            if size < maxBitsSize:
                return ovlHashBitsDict[size]
