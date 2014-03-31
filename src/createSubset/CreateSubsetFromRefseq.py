import sys
from model import Pool, Command
from assembly import WgsAssembler

class Main(object):
    
    def main(self):
        pool = self.parseArguments()
        self.createSubset(pool)
        WgsAssembler.WgsAssembler().doAssembly(pool)
        
    def parseArguments(self):
        if len(sys.argv) == 1:
            self.printHelp()
        elif len(sys.argv)==2:
            return Pool.Pool(sys.argv[1])
        else:
            self.printHelp()
            
    def printHelp(self):
        print("Usage: CreateSubset.py <configuration file>")
        sys.exit()
        
    def createSubset(self, pool):
        for lib in pool.libs:
            BwaCommand(lib)#.execute()
            SamToFastqCommand(lib)#.execute()
    
class BwaCommand(Command.Command):
    def setCommand(self):
        self.outBam = self.collection.outputDir+self.collection.libName+".bam"
        
        self.addArg("bwa mem  24 -U 2")
        self.addArg("-t " + self.collection.config.getGlobalOption("maxThreads"))
        self.addArg("/local/work/jetse/Parasponia/data/fastAssembly/paraLongestSequences.fasta")
        self.addArg(self.collection.forward)
        self.addArg(self.collection.reversed)
        self.addArg("| samtools view")
        self.addArg("-b -S -F 12 -h -")
        self.addArg(">")
        self.addArg(self.outBam)
    
    def updateStatus(self):
        self.collection.bam = self.outBam

class SamToFastqCommand(Command.Command):
    def setCommand(self):
        self.inBam = self.collection.bam
        self.outForward = self.collection.outputDir+self.collection.libName+"_1.fastq"
        self.outReversed = self.collection.outputDir+self.collection.libName+"_2.fastq"
        
        self.addArg("java -jar /local/work/jetse/programs/picard-tools-1.107/SamToFastq.jar")
        self.addArg("INPUT=" + self.inBam)
        self.addArg("FASTQ=" + self.outForward)
        self.addArg("SECOND_END_FASTQ=" + self.outReversed)
          
    def updateStatus(self):
        self.collection.forward = self.outForward
        self.collection.reversed = self.outReversed

if __name__ == '__main__':
    Main().main()