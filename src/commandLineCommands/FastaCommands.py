import Command, os
from configuration import Configuration

class CegmaCommand(Command.Command):
    
    def setCommand(self):
        name = "cegma"
        self.addArg("cegma")
        self.addArg("-g " + self.genome)
        if self.coreGenes !=None:
            self.addArg("-p " + self.coreGenes)
            name = os.path.basename(os.path.splitext(self.coreGenes)[0])
        self.addArg("-o " + self.outputDir + name)
        self.addArg("-T " + Configuration.instance.getGlobalOption("maxThreads"))
        
        self.outputFile = self.outputDir + name + ".completeness_report"
        
class BlastCommand(Command.Command):
    
    def setCommand(self):
        self.outputFile = self.outputDir + "contaminationCheck.txt"
        self.addArg("blastn")
        self.addArg("-db " + self.db)
        self.addArg("-query " + self.fastaFile)
        self.addArg("-out " + self.outputFile)
        self.addArg("-outfmt \"6 qseqid sseqid evalue bitscore sgi sacc staxids stitle\"")
        self.addArg("-num_alignments 2")
        
class GenBlastA(Command.Command):
    
    def setCommand(self):
        self.outputFile = self.outputDir + os.path.basename(os.path.splitext(self.proteins)[0]) + ".gbo"
        self.report = self.referenceGenome + ".gba.report"
        self.blastOutput = self.referenceGenome + ".blast"
        self.addArg("genBlastA")
        self.addArg("-q " + self.proteins)
        self.addArg("-t " + self.referenceGenome)
        self.addArg("-o " + self.outputFile)