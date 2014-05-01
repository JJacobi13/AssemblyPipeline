from configuration import Configuration
from commandLineCommands import Command

class JellyFishCount(Command.Command):
    """
    The class JellyFishCounts executes jellyfish count to create a file genSize_0 with count information for each different kmer.
    """
    def checkInput(self):
        for lib in self.pool.libs:
            if hasattr(lib, "forward") == False:
                raise ValueError("Library has no forward reads, which are required for executing the jellyfish command")
    
    def setCommand(self):
        self.outputFile = self.outputDir + "genSize_0"
        self.addArg("jellyfish count")
        self.addArg("-m " + Configuration.instance.getGlobalOption("kmer"))
        self.addArg("-s " + str(int(Configuration.instance.getGlobalOption("maxMem")) * 1000))
        self.addArg("-t " + Configuration.instance.getGlobalOption("maxThreads"))
        self.addArg("-o " + self.outputDir + "genSize")
        self.addArg("-C")
        for lib in self.pool.libs:
            self.addArg(lib.forward)
            if lib.reversed != None:
                self.addArg(lib.forward)
                
class JellyFishDump(Command.Command):
    """
    The class JellyFishDump converts the kmer information to a human readable format.
    """
    def setCommand(self):
        self.outputFile = self.outputDir + "kmer.counts"
        self.addArg("jellyfish dump " + self.inputFile + " > " + self.outputFile)
    
class JellyFishStats(Command.Command):
    """
    This class calculates some statistics about the the kmers with jellyfish stats
    """
    def checkInput(self):
        if hasattr(self, "jellyFishCountsFile") == False:
            raise ValueError("The JellyFishStats command requires a JellyfishCounts file")
    
    def setCommand(self):
        self.outputFile = self.outputDir + "genomeSize_stats.txt"
        self.addArg("jellyfish stats")
        self.addArg(" -o " + self.outputFile)
        self.addArg(self.jellyFishCountsFile)
#             self.addArg(self.outputDir + "genSize_0")
        
class JellyFishHisto(Command.Command):  
    """
    This class generates a histogram of all kmer counts.
    """
    def checkInput(self):
        if hasattr(self, "jellyFishCountsFile") == False:
            raise ValueError("The JellyFishHisto command requires a jellyFishCountsFile file")
        
    def setCommand(self):
        self.outputFile = self.outputDir + "genomeSize_histo.histo"
        self.addArg("jellyfish histo")
        self.addArg(" -h " + str((int(Configuration.instance.getGlobalOption("expCoverage")) * 2)))
        self.addArg(" -t " + Configuration.instance.getGlobalOption("maxThreads"))
        self.addArg(" -o " + self.outputFile)
        self.addArg(self.jellyFishCountsFile)            