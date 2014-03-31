import Command

class MergeCommand(Command.Command):
    
    def setCommand(self):
        self.outputFile = self.outputDir + "mergedFiles.fastq"
        self.addArg("cat")
        for fastqFile in self.fastqFiles:
            self.addArg(fastqFile)
        self.addArg("> " + self.outputFile)
        