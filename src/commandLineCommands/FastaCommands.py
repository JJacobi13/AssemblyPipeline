import Command
from configuration import Configuration

class CegmaCommand(Command.Command):
    
    def setCommand(self):
        self.outputFile = self.outputDir + "cegma.completeness_report"
        self.addArg("cegma")
        self.addArg("-g " + self.genome)
        self.addArg("-o " + self.outputDir + "cegma")
        self.addArg("-T " + Configuration.instance.getGlobalOption("maxThreads"))
        
class BlastCommand(Command.Command):
    
    def setCommand(self):
        self.outputFile = self.outputDir + "contaminationCheck.txt"
        self.addArg("blastn")
        self.addArg("-db " + self.db)
        self.addArg("-query " + self.fastaFile)
        self.addArg("-out " + self.outputFile)
        self.addArg("-outfmt \"6 qseqid sseqid evalue bitscore sgi sacc staxids stitle\"")
        self.addArg("-num_alignments 2")
        
        