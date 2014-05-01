from commandLineCommands import Command, Mappers
from configuration import Configuration
import os

class ContaminationFiltering(object):
    def filterContamination(self, outputDir, forwardReads, reversedReads, insertSize, contFasta):
#         phixGenome = os.path.dirname(__file__) + "/phixDb/phix.fasta"
        contamination = os.path.splitext(os.path.basename(contFasta))[0]
        filtered = BowtieFilter(outputDir, refGenome = contFasta, forward = forwardReads, reversed = reversedReads, insertSize=insertSize, contamination= contamination).execute()
        if reversedReads !=None:
            return [filtered, filtered.replace(".1.fastq",".2.fastq")]
        return [filtered, None]
    
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
        self.baseoutputFile = self.outputDir+self.contamination + "Filtered.fastq"
        
        self.addArg("bowtie2")
        self.addArg("-p " + Configuration.instance.getGlobalOption("maxThreads"))
        self.addArg("-x " + self.bowtieIndex)
        if self.reversed == None:
            self.addArg("--un " + self.baseoutputFile)
            self.addArg("-U " + self.forward)
            self.outputFile = self.baseoutputFile
        else:
            self.reversedFile = self.outputDir+self.contamination + "Filtered.2.fastq"
            self.outputFile = self.outputDir+self.contamination + "Filtered.1.fastq"
            self.addArg("--un-conc " + self.baseoutputFile)
            self.addArg("-I " + str(int(int(self.insertSize) * 0.5)))
            self.addArg("-X " + str(int(self.insertSize)*2))
            self.addArg("-1 " + self.forward)
            self.addArg("-2 " + self.reversed)
            
        self.addArg("-S /dev/null")