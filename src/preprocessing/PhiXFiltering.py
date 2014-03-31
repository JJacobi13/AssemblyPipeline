from model import Command
import os

class PhiXFiltering(object):
    pass

class Bowtie(Command.Command):
    def checkInput(self):
        if os.path.isfile(self.collection.phixDb + ".1.bt2") == False:
            BowtieIndex().execute()
            
    def setCommand(self):
        self.outputFile = self.collection.outputDir + self.collection.libName
        
        self.addArg("bowtie2")
        self.addArg("-p " + self.collection.config.getGlobalOption("maxThreads"))
        self.addArg("-x " + self.collection.bowtieIndex)
        if hasattr(self.collection, "reversed"):
            self.addArg("-1 " + self.collection.forward)
            self.addArg("-2 " + self.collection.reversed)
        else:
            self.addArg("-U " + self.collection.forward)
        self.addArg("-S " + self.outputFile)
    
    def updateStatus(self):
        self.collection.samFile = self.outputFile

class BowtieIndex(Command.Command):
    
    def checkInput(self):
        pass
            
    def setCommand(self):
        self.bowtieIndex = os.path.dirname(__file__) + "/phixDb/phix"
        self.outputFile = self.bowtieIndex + ".1.bt2"
        
        self.addArg("bowtie2-build")
        self.addArg(os.path.dirname(__file__) + "/phixDb/phix.fasta")
        self.addArg(self.bowtieIndex)
    
    def updateStatus(self):
        self.collection.phixDb = self.bowtieIndex
    