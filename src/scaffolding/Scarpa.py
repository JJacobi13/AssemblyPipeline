from commandLineCommands import Command
from utils import FastqUtils
from configuration import Configuration
import os

class Scarpa():
    
    def doScaffolding(self, outputDir, assembly, pool):
        for lib in pool.libs:
            lib.mergedFile = outputDir + lib.libName + "Merged.fastq"
            if os.path.exists(lib.mergedFile) == False or Configuration.instance.getGlobalOption("overwrite") == 1:
                FastqUtils.mergeFastqFiles(lib.forward, lib.reversed, lib.mergedFile)
        scarpaProcesser = ScarpaProcess(outputDir, contigs=assembly, libs=pool.libs)
        scarpaContigs = scarpaProcesser.execute()
        libInfo = scarpaProcesser.libInfo
        mergedFiles = []
        for lib in pool.libs:
            mergedFiles.append(lib.mergedFile)
        samFile = BowtieOne(outputDir, contigs=scarpaContigs, mergedFiles=mergedFiles).execute()
        mapFile = ScarpaParser(outputDir, samFile=samFile).execute()
        self.scaffolds = ScarpaCommand(outputDir, libInfo = libInfo, contigs = scarpaContigs, mappedFile=mapFile).execute()
        
        return self.scaffolds        
        
class ScarpaProcess(Command.Command):

    def setCommand(self):
        self.outputFile = self.contigs + ".scarpa.fa"
        self.libInfo = self.contigs + ".scarpa.info"
        
        self.addArg("scarpa_process")
        for lib in self.libs:
            self.addArg("-f " + lib.mergedFile)
            self.addArg("-i " + lib.insertSize)
            lib.mergedFile = lib.mergedFile + ".scarpa.fq"
        self.addArg("-c " + self.contigs)
        
class BowtieOne(Command.Command):
    
    def checkInput(self):
        bowtieIndexer = BowtieOneIndexer(self.outputDir, contigs=self.contigs)
        bowtieIndexer.execute()
        self.bowtieIndex = bowtieIndexer.index
    
    def setCommand(self):
        self.outputFile = os.path.splitext(self.contigs)[0] + ".sam"
        self.addArg("bowtie -p " + Configuration.instance.getGlobalOption("maxThreads"))
        self.addArg("--sam")
        if FastqUtils.determineQuality(self.mergedFiles[0]) == 64:
            self.addArg("--phred64-quals")
        self.addArg(self.bowtieIndex)
        self.addArg(",".join(self.mergedFiles))
        self.addArg(self.outputFile)
        
class ScarpaParser(Command.Command):
    def setCommand(self):
        self.outputFile = os.path.splitext(self.samFile)[0] + ".map"
        self.addArg("cat " + self.samFile + " |")
        self.addArg("scarpa_parser > " + self.outputFile)
        
class BowtieOneIndexer(Command.Command):
    def setCommand(self):
        self.index = os.path.splitext(self.contigs)[0]
        self.outputFile = self.index + ".1.ebwt"
        self.addArg("bowtie-build " + self.contigs)
        self.addArg(self.index)
        
class ScarpaCommand(Command.Command):
    def setCommand(self):
        self.outputFile = self.outputDir + "scaffolds.fasta"
        self.addArg("scarpa -c " + self.contigs)
        self.addArg("-l " + self.libInfo)
        self.addArg("-i " + self.mappedFile)
        self.addArg("-o " + self.outputFile)
        
if __name__ == '__main__':
    from model import Pool
    import sys, logging
    
    logging.basicConfig(filename="pipeline.log", format="%(asctime)-25s%(levelname)-12s%(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p")
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter("%(asctime)-25s%(message)s", "%I:%M:%S %p"))
    logging.getLogger().addHandler(console)

    Configuration.instance.parseIni(sys.argv[1])
    pool = Pool.Pool(Configuration.instance.getGlobalOption("outputDirectory"))
    pool.createLibs(Configuration.instance.libNames)
    
    Scarpa().doScaffolding(pool.outputDir, sys.argv[2], pool)
    
    