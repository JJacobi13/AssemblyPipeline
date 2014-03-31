from commandLineCommands import Command
import os

class SffToFastqConverter(Command.Command):
    
    def setCommand(self):
        self.outputFile = self.outputDir + os.path.basename(os.path.splitext(self.sffFile)[0]) + ".fastq"
        self.addArg("sff2fastq")
        self.addArg("-o " + self.outputFile)
        self.addArg(self.sffFile)
    