import Command
from configuration import Configuration
import os

class Bowtie(Command.Command):
    """
    The bowtie command to map DNA reads against a reference genome with bowtie
    """
    def checkInput(self):
        if hasattr(self,"bowtieIndex") == False:
            indexer = BowtieIndex(self.outputDir, refGenome=self.refGenome)
            indexer.execute()
            self.bowtieIndex = indexer.bowtieIndex
            
    def setCommand(self):
        self.outputFile = self.outputDir+"mapped.sam"
        
        self.addArg("bowtie2")
        self.addArg("-p " + Configuration.instance.getGlobalOption("maxThreads"))
        self.addArg("-x " + self.bowtieIndex)
        self.addArg("")
        if hasattr(self, "fastqFile"):
            self.addArg("-U " + self.fastqFile)
        else:
            self.addArg("-I " + str(int(int(self.insertSize) * 0.5)))
            self.addArg("-X " + str(int(self.insertSize)*2))
            self.addArg("-1 " + self.forward)
            self.addArg("-2 " + self.reversed)
        self.addArg("-S " + self.outputFile)

class BowtieIndex(Command.Command):
    """
    The class BowtieIndex creates a bowtie index for a reference genome (or assembly)
    """                
    def setCommand(self):
        self.bowtieIndex = os.path.splitext(self.refGenome)[0] 
        self.outputFile = self.bowtieIndex + ".1.bt2"
        
        self.addArg("bowtie2-build")
        self.addArg(self.refGenome)
        self.addArg(self.bowtieIndex)

    
class Tophat(Command.Command):
    """
    The Tophat command executes tophat to map RNA reads against a reference genome.
    """
    def checkInput(self):
        if hasattr(self,"bowtieIndex") == False:
            indexer = BowtieIndex(self.outputDir, refGenome=self.refGenome)
            indexer.execute()
            self.bowtieIndex = indexer.bowtieIndex
    
    def setCommand(self):
        self.outputFile = self.outputDir + os.path.basename(os.path.splitext(self.forwardRna)[0]) + "/accepted_hits.bam"
        self.addArg("tophat")
        self.addArg("-o " + os.path.dirname(self.outputFile))
        self.addArg(self.bowtieIndex)
        self.addArg(self.forwardRna)
        self.addArg(self.reversedRna)
        
class BwaSampe(Command.Command):
    
    def checkInput(self):
        BwaIndex(self.outputDir, fastaFile=self.refGenome).execute()
        
    def setCommand(self):
        self.outputFile = self.outputDir + os.path.basename(os.path.splitext(self.forwardFastq)[0]) + ".sam"
        self.addArg("bwa sampe " + self.refGenome)
        self.addArg(BwaAln(self.outputDir, fastqFile = self.forwardFastq, refGenome=self.refGenome).execute())
        self.addArg(BwaAln(self.outputDir, fastqFile = self.reversedFastq, refGenome=self.refGenome).execute())
        self.addArg(self.forwardFastq)
        self.addArg(self.reversedFastq)
        self.addArg(self.outputFile)
        
class BwaAln(Command.Command):
    
    def checkInput(self):
        BwaIndex(self.outputDir, fastaFile=self.refGenome).execute()
        
    def setCommand(self):
        self.outputFile = self.outputDir + os.path.basename(os.path.splitext(self.fastqFile)[0]) + ".sai"
        self.addArg("bwa aln ")
        self.addArg("-t " + Configuration.instance.getGlobalOption("maxThreads"))
        self.addArg(self.refGenome)
        self.addArg(self.fastqFile)
        self.addArg(self.outputFile)
        
class BwaIndex(Command.Command):
    
    def setCommand(self):
        self.outputFile = self.fastaFile+".pac"
        self.addArg("bwa index " + self.fastaFile)
        
    