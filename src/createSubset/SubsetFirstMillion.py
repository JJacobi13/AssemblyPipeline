from commandLineCommands import Command
import os

class SubsetFirstMillion(Command.Command):
    def getDescription(self):
        return "Create a subset of first million reads of a  fastq file"
    
    def setCommand(self):
        self.outputFile = self.outputDir + "sub_" + os.path.basename(self.fastqFile)
        self.addArg("head")
        self.addArg("-n 4000000")
        self.addArg(self.fastqFile)
        self.addArg("> " + self.outputFile)
    
    