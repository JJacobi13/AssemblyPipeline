from commandLineCommands import Command, Mappers
from configuration import Configuration
import os, logging

class PhiXFiltering(object):
    def filterPhix(self, outputDir, forwardReads, reversedReads, insertSize):
        phixGenome = os.path.dirname(__file__) + "/phixDb/phix.fasta"
        filtered = BowtieFilter(outputDir, refGenome = phixGenome, forward = forwardReads, reversed = reversedReads, insertSize=insertSize, contamination= "PhiX").execute()
        splittedExt = os.path.splitext(filtered)
        if reversedReads !=None:
            return [splittedExt[0] + ".1" + splittedExt[1], splittedExt[0] + ".2" + splittedExt[1]]
        return filtered
    
class BowtieFilter(Command.Command):
    """
    The bowtie command to map DNA reads against a reference genome with bowtie
    """
    
    def getDescription(self):
        return "Removed contamination of " + self.contamination + " with bowtie version 2.1.0 (WARNING: hardcoded)"
    
    def checkInput(self):
        indexer = Mappers.BowtieIndex(self.outputDir, refGenome=self.refGenome)
        indexer.execute()
        self.bowtieIndex = indexer.bowtieIndex
            
    def setCommand(self):
        self.baseoutputFile = self.outputDir+"phixFiltered.fastq"
        
        self.addArg("bowtie2")
        self.addArg("-p " + Configuration.instance.getGlobalOption("maxThreads"))
        self.addArg("-x " + self.bowtieIndex)
        if self.reversed == None:
            self.addArg("--un " + self.baseoutputFile)
            self.addArg("-U " + self.forward)
            self.outputFile = self.baseoutputFile
        else:
            self.reversedFile = self.outputDir+"phixFiltered.2.fastq"
            self.outputFile = self.outputDir+"phixFiltered.1.fastq"
            self.addArg("--un-conc " + self.baseoutputFile)
            self.addArg("-I " + str(int(int(self.insertSize) * 0.5)))
            self.addArg("-X " + str(int(self.insertSize)*2))
            self.addArg("-1 " + self.forward)
            self.addArg("-2 " + self.reversed)
            
        self.addArg("-S /dev/null")