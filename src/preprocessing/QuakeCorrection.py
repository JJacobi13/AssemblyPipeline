from commandLineCommands import JellyFish, Command
from configuration import Configuration

from model import Pool, Library
import logging

class QuakeCorrection():
    
    def doCorrection(self, pool, cutoff):
        libs = []
        for lib in pool.libs:
            if lib.readlen < 125:
                libs.append(lib)
            else:
                logging.info("Skipping quake correction for library " +lib.libName+", quake only works with reads shorter than 125bp") 

        countsFile = QmerCounter(pool.outputDir, libs=libs).execute()
        QuakeCorrect(lib.outputDir, libs=libs, jellyfishFile = countsFile, cutoff = cutoff).execute()
        
        
class QuakeCorrect(Command.Command):
    
    def setCommand(self):
        self.outputFile = "aFile"
        self.addArg("correct -f " + self.getFastqFileNames())
        self.addArg("-k " + Configuration.instance.getGlobalOption("kmer"))
        self.addArg("-c " + str(self.cutoff))
        self.addArg("-m " + self.jellyfishFile)
        self.addArg("-p " + Configuration.instance.getGlobalOption("maxThreads"))
        
    def getFastqFileNames(self):
        kmerFile = self.outputDir + "kmerFiles.txt"
        with open(kmerFile, "w") as writer:
            for lib in self.libs:
                writer.write(lib.forward)
                if lib.reversed !=None:
                    writer.write(" " + lib.reversed + "\n")
        return kmerFile
                    
        
class QmerCounter(Command.Command):
    def setCommand(self):
        self.outputFile = self.outputDir + "kmer.counts"
        self.addArg("cat ")
        for lib in self.libs:
            self.addArg(lib.forward)
            if lib.reversed != None:
                self.addArg(lib.forward)
        self.addArg("|")
        self.addArg("count-qmers -k " + Configuration.instance.getGlobalOption("kmer"))
        self.addArg("> " + self.outputFile)
        
        
pool = Pool.Pool(".")
lib = Library.Library(pool,"test")
lib.forward = "gist_pf_trimmed_1.fastq"
lib.reversed = "gist_pf_trimmed_2.fastq"
pool.libs = [lib]
Configuration.instance.setOption("maxThreads", "20")
Configuration.instance.setOption("overwrite", "0")
Configuration.instance.setOption("kmer", "17")
Configuration.instance.setOption("maxMem", "200000000")
QuakeCorrection().doCorrection(pool, 32)
