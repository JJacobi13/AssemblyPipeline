import Command, os

class BamMerger(Command.Command):
    
    def setCommand(self):
        self.outputFile=self.outputDir + "mergedBams.bam"
        self.addArg("samtools merge -n")
        self.addArg(self.outputFile)
        for bamFile in self.bamFiles:
            self.addArg(bamFile)        

class SamToBamConverter(Command.Command):
    
    def setCommand(self):
        self.outputFile = os.path.splitext(self.samFile)[0] + ".bam"
        self.addArg("samtools view -bS " + self.samFile + "|")
        self.addArg("samtools sort - " + os.path.splitext(self.outputFile)[0])

class SamtoolsMpileup(Command.Command):
    
    def checkInput(self):
        if hasattr(self, "bamFile") == False:
            self.bamFile = SamToBamConverter(self.samFile).execute()
    
    def setCommand(self):
        self.outputFile = os.path.splitext(self.bamFile)[0] + ".vcf"
        self.addArg("samtools mpileup -Sguf " + self.fastaFile + " " +self.bamFile+" |")
        self.addArg("bcftools view -vcg - > " + self.outputFile)