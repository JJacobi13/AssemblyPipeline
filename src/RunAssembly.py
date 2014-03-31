import sys, os, copy
from model import Pool
from genomeSizeEstimation import GenomeSizeEstimation
from qualityControl import AssemblyControl,FastqControl, InsertSizeChecker, BlastScanner
from qualityControl.Reporter import Reporter
from configuration import Configuration
from preprocessing import FastqMcfTrimming, SffToFastqConverter
from assembly import WgsAssembler
from commandLineCommands import FastaCommands, FastqCommands
from scaffolding import PBJelly
import logging
from utils import DirUtils
# from qualityControl import FastqSmallReport

class RunAssembly(object):
    """
    The RunAssembly class is the main flow of the full pipeline. This class makes the decisions which programs 
    will be executed.
    """
    
    def __init__(self):
        """
        The constructor of the RunAssembly class makes the highest level decisions of the pipeline.
        """
        pool = self.parseArguments()
        for lib in pool.libs:
             
            lib.forward = DirUtils.fileRegexToList(lib.forward)
            if lib.format == "sff":
                for idx, sffFile in enumerate(lib.forward):
                    lib.forward[idx] = SffToFastqConverter.SffToFastqConverter(lib.outputDir, sffFile=sffFile).execute()
            logging.info("Creating small fastq report of " + lib.libName)
#             smallReportGenerator = FastqSmallReport.FastqSmallReport()
#             smallReportGenerator.createSmallReport(lib.forward,lib.libName)
#             Reporter.instance.objects.append(smallReportGenerator)
             
            lib.forward = FastqCommands.MergeCommand(lib.outputDir, fastqFiles=lib.forward).execute()
             
        self.doGenomeSizeEstimation(pool)    
        
        for lib in pool.libs:
            if lib.sequencingPlatform == "illumina":
                self.illuminaPreprocess(lib)
            elif lib.sequencingPlatform == "454":
                reports = self.createFastqReport(lib, lib.outputDir + "raw_qc/")
                lib.forward = FastqMcfTrimming.FastqTrimmer(lib.outputDir, fastqReports=reports, forward=lib.forward,noTrim=True).execute()
                reports = self.createFastqReport(lib, lib.outputDir + "preprocessed/")
        self.doAssembly(pool)   
        for lib in pool.libs:
            if lib.sequencingPlatform == "pacbio":
                self.assembly = PBJelly.PBJelly(pool.outputDir + "scaffolds/", assembly=self.assembly, reads=lib.forward).execute()
        self.contaminationCheck(self.assembly, pool.outputDir)
        
        logging.info("Creating assembly stats")
        assemblyController = AssemblyControl.AssemblyStatistics()
        assemblyController.AssemblyStatisticsOfPipeline(pool.outputDir + "statistics/", pool, self.assembly)
        
        Reporter.instance.createReport(pool.outputDir + "report/")

        
    def parseArguments(self):
        """
        The method parseArguments parses the given arguments. Now only extracts the first argument as configuration file.
        #TODO: parse all other key-value arguments which also can be given in the configuration file.
        """
        if len(sys.argv) == 1:
            self.printHelp()
        elif len(sys.argv)==2:
            Configuration.instance.parseIni(sys.argv[1])
            pool = Pool.Pool(Configuration.instance.getGlobalOption("outputDirectory"))
            pool.createLibs(Configuration.instance.libNames)
            return pool
        else:
            self.printHelp()
            
    def printHelp(self):
        """
        This method prints the help, very simple now...
        """
        print("Usage: python RunAssembly.py <configuration file>")
        sys.exit()
        
    def contaminationCheck(self, assembly, outputDir):
        """
        This method creates objects to execute the contamination check of a given assembly. 
        The output is written to the given output directory.
        """
        blastFile = FastaCommands.BlastCommand(outputDir, fastaFile=self.assembly, db=Configuration.instance.getGlobalOption("nrDb")).execute()
        blastScanner = BlastScanner.BlastScanner()
        blastScanner.scan_results(blastFile)
        return [blastScanner]
        
    def doAssembly(self, pool):
        """
        The method doAssembly creates all objects to execute a wgs assembly. Afterwards the insert sizes of all pe and
        mp libraries are estimated.
        """
        logging.info("Executing assembly")
        assembler = WgsAssembler.WgsAssembler()
        self.assembly = assembler.doAssembly(pool.outputDir + "assembly/", pool)
          
        for lib in pool.libs:
            if lib.reversed == None:
                continue
            logging.info("Calculating insert sizes for " + lib.libName)
            insertSizeChecker = InsertSizeChecker.InsertSizeChecker()
            insertSizeChecker.checkInsertSize(lib.outputDir, lib.rawForward, lib.rawReversed, self.assembly, lib.libName, lib.insertSize)
        
    def doGenomeSizeEstimation(self, pool):
        """
        This method creates the objects to do a genome size estimation.
        """
        genomeSizeEstimator = GenomeSizeEstimation.GenomeSizeEstimation()
        genomeSizeEstimator.doGenomeSizeEstimation(pool.outputDir + "genomeSizeEstimation/", pool)
        
    def illuminaPreprocess(self, lib):
        """
        This method manages the illumina preprocessing by creating fastqMcfTrimming objects to do the trimming/adapter
        removal/duplicate removal.
        """
        fastqReports = self.createFastqReport(lib, lib.outputDir + "raw_qc/")
        #Run fastqmcf to trim the reads
        logging.info("Trimming " + lib.libName)
        if lib.rawReversed !=None:
            fastqTrimmer = FastqMcfTrimming.FastqTrimmer(lib.outputDir, fastqReports=fastqReports, forward=lib.rawForward,reversed=lib.rawReversed)
        else:
            fastqTrimmer = FastqMcfTrimming.FastqTrimmer(lib.outputDir, fastqReports=fastqReports, forward=lib.rawForward)
        lib.forward = fastqTrimmer.execute()
        if lib.rawReversed !=None:
            lib.reversed = fastqTrimmer.outReversed
           
        #Create fastq report of the preprocessed data
        self.createFastqReport(lib, lib.outputDir + "preprocessed_qc/")
        
    
    def createFastqReport(self, lib, outDir):
        """
        This method creates the fastqc report object which executes fastqc.
        """
        logging.info("Creating fastq report for " + lib.libName)
        fastqController = FastqControl.FastqReportCreator(outDir, forward=lib.rawForward, status="forward raw reads of " + lib.libName)
        forwardQcIndex = fastqController.execute()
        if lib.rawReversed !=None:
            rawFastqControllerRev = FastqControl.FastqReportCreator(outDir, forward=lib.rawReversed, status="reversed raw reads of " + lib.libName)
            reversedQcIndex = rawFastqControllerRev.execute()
            return[forwardQcIndex,reversedQcIndex]
        return[forwardQcIndex]
        
if __name__ == '__main__':
    RunAssembly()