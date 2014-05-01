import Command

class MergeCommand(Command.Command):
    """
    The mergecommand merges all given fastqfiles with cat into a single file.
    """
    
    def getDescription(self):
        return "Concatenation of all fastq files"
    
    def setCommand(self):
        self.outputFile = self.outputDir + self.direction +"Merged.fastq"
        self.addArg("cat")
        for fastqFile in self.fastqFiles:
            self.addArg(fastqFile)
        self.addArg("> " + self.outputFile)
        