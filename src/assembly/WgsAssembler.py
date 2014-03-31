import sys, os, math, logging
from model import Pool
from commandLineCommands import Command
from configuration import Configuration
'''
Created on Feb 14, 2014

@author: Jetse
'''
class WgsAssembler(object):
    """
    The WgsAssembler class creates all commands to run the wgs assembler.
    """
    
    def main(self):
        """
        The main method of the WgsAssembler has to be called when only the wgs assembler has to be executed. The configuration file has to be given as the only argument of this script
        """
        pool = self.parseArguments()
        self.doAssembly(pool)
        
    def doAssembly(self, outputDir, pool):
        """
        The method doAssembly creates the commands to run the wgs assembler
        :param pool: The pool of libraries to execute the wgs assembler on
        :type pool: a Pool object
        """
        frgFiles = []
        for lib in pool.libs:
            if lib.sequencingPlatform != "pacbio":
                frgFiles.append(FastqToCaCommand(outputDir,libName=lib.libName, insertSize=lib.insertSize,stdev=lib.stdev,readlen=lib.readlen,type=lib.type,forward=lib.rawForward,reversed=lib.rawReversed, sequencingPlatform=lib.sequencingPlatform).execute())
        assembly = RunCaCommand(outputDir, frgFiles=frgFiles).execute()
        return assembly
        
    def parseArguments(self):
        """
        The method parseArguments creates a pool object with the first given parameter to the script. This parameter is expected to be a configuration file.
        """
        if len(sys.argv) == 1:
            self.printHelp()
        elif len(sys.argv)==2:
            pool = Pool.Pool(sys.argv[1])
            return pool
        else:
            self.printHelp()
            
    def printHelp(self):
        """
        The method printHelp print the help and exits the program
        """
        print("Usage: WgsAssembler <configuration file>")
        sys.exit()
        
class SpecFile(object):
    """
    The specfile object creates a specfile for the wgs assembler based on the number of threads and threads the user wants to use
    """
    def __init__(self, outputDir):
        self.fileName = outputDir + "wgs.spec"
        maxMem = int(Configuration.instance.getGlobalOption("maxMem"))
        maxThr = int(Configuration.instance.getGlobalOption("maxThreads"))
        with open(self.fileName, "w") as specWriter:
            specWriter.write("#Default allowed error rates for Illumina\n")
            specWriter.write("utgErrorRate = 0.03\n")
            specWriter.write("utgErrorLimit = 6.5\n")
            specWriter.write("cnsErrorRate = 0.06\n")
            specWriter.write("cgwErrorRate = 0.10\n")
            specWriter.write("ovlErrorRate = 0.06\n")
            specWriter.write("merSize=22\n")
            specWriter.write("unitigger=bogart\n")
            
            specWriter.write("#Don't use grid\n")
            specWriter.write("useGrid = 0\n")
            specWriter.write("scriptOnGrid = 0\n")
            specWriter.write("frgCorrOnGrid = 0\n")
            specWriter.write("ovlCorrOnGrid = 0\n")
            
            specWriter.write("#Memory settings for the given maximum memory/threads\n")
            specWriter.write("merylMemory = " + str(int(maxMem/1000)) + "\n") #merylMemory is in mb
            specWriter.write("merylThreads = " + str(maxThr) + "\n")
            
            specWriter.write("ovlStoreMemory = " + str(int(maxMem/1000)) + "\n") #ovlStoreMemory is in mb
            
            hashBits = self.getOvlHashBits(maxMem*0.66, maxThr)
            specWriter.write("ovlThreads = " + str(self.ovlThreads) + "\n")
            specWriter.write("ovlConcurrency = " + str(self.ovlConcurrency) + "\n")
            specWriter.write("ovlHashBits = " + str(hashBits) + "\n")

            specWriter.write("ovlHashBlockLength = " + str(maxMem*0.33/self.ovlConcurrency*100) + "\n")
            specWriter.write("ovlRefBlockSize = 20000000\n")
            
            specWriter.write("ovlCorrBatchSize = 1000000\n")
            specWriter.write("ovlCorrConcurrency = "+ str(int(math.floor(min(maxThr,int(math.floor(maxMem/1000))))))+"\n")
            
            frgCorrThreads = self.ovlThreads
            specWriter.write("frgCorrThreads = " + str(frgCorrThreads) + "\n")
            frgCorrConcurrency = self.ovlConcurrency
            specWriter.write("frgCorrConcurrency = " + str(frgCorrConcurrency) + "\n")
            frgCorrBatchSize = (100/float(1300))*maxMem
            specWriter.write("frgCorrBatchSize = " + str(int(frgCorrBatchSize)) + "\n")
            
            specWriter.write("merOverlapperThreads = " + str(maxThr) + "\n")
            #edit merOverlapperSeedBatchSize when documentation available
            specWriter.write("cnsConcurrency = " + str(maxThr) + "\n")
            
            specWriter.write("mbtThreads = " + str(frgCorrThreads) + "\n")
            specWriter.write("mbtConcurrency = " + str(frgCorrConcurrency) + "\n")
            #relationship mbtBatchSize and memory/speed?
            
            
    def getOvlHashBits(self, maxMem,maxThr):
        self.ovlThreads = int(math.floor(math.sqrt(maxThr)))
        self.ovlConcurrency = int(math.floor(maxThr/self.ovlThreads))
        ovlHashBitsDict = {54:18,108:19,216:20,432:21,864:22,1728:23,3456:24,6912:25,13824:26,27648:27,55296:28,110592:29,221184:30}
        maxBitsSize = maxMem/self.ovlConcurrency/1000
        for size in sorted(ovlHashBitsDict, reverse=True):
            if size < maxBitsSize:
                return ovlHashBitsDict[size]
        
class RunCaCommand(Command.Command):
    """
    The RunCaCommand runs the wgs assembler
    """  
    def setCommand(self):
        self.assemblyPrefix = "assembly"
        self.outputFile = self.outputDir + "9-terminator/"+self.assemblyPrefix+".ctg.fasta"
        if not os.path.isdir(self.outputDir):
            os.makedirs(self.outputDir)
        specFile = SpecFile(self.outputDir)
        
        self.addArg("runCA")
        self.addArg("-d " + self.outputDir)
        self.addArg("-p " + self.assemblyPrefix)
        self.addArg("-s " + specFile.fileName)
        for frgFile in self.frgFiles:
            self.addArg(frgFile)
        
class FastqToCaCommand(Command.Command):
    """
    The FastqToCaCommand prepares the input for the wgs assembler with the command fastqToCA
    """

    def setCommand(self):
        self.outputFile = self.outputDir+self.libName + ".frg"    
            
        self.addArg("fastqToCA")
        if self.reversed != None:
            self.addArg("-insertsize " + self.insertSize + " " + self.stdev)
            self.addArg(self.getOrientation())
            self.addArg("-mates "+self.forward+","+self.reversed)
        else:
            self.addArg("-reads " + self.forward)
        self.addArg("-libraryname " + self.libName)
        self.addArg("-technology " + self.getTech())
        
        
        self.addArg(">")
        self.addArg(self.outputFile)

    def getOrientation(self):
        """
        The method getOrientation returns the orientation of the reads.
        TODO: implement orientation for other types of data
        """
        logging.warn("Orientation is determined for Illumina reads: mp reads are outwards (forward: 3'->5', reversed: 5'->3'), PE reads are inwards (forward: 5'->3', reversed: 3'->5')")
        if self.type == "mp":
            return "-outtie"
        elif self.type == "pe":
            return "-innie"
    
    def getTech(self):
        """
        The method getType returns the type of the reads.
        TODO: implement orientation for other types of data
        """
        if self.readlen > 160:
            return "illumina-long"
        elif self.sequencingPlatform == "illumina":
            return "illumina"
        elif self.sequencingPlatform == "454":
            return "454"

        
"""
If this script is called from the command-line, execute the main...
"""
if __name__ == '__main__':
    WgsAssembler().main()