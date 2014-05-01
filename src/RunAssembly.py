import sys, os
from model import Pool
from genomeSizeEstimation import GenomeSizeEstimation
from qualityControl import AssemblyControl,FastqControl, InsertSizeChecker, BlastScanner, FastqSmallReport
from qualityControl.Reporter import Reporter, LaTeX
from configuration import Configuration
from preprocessing import FastqMcfTrimming, SffToFastqConverter, ContaminationFiltering, QuakeCorrection
from assembly import WgsAssembler, AllpathsAssembler
from commandLineCommands import FastaCommands, FastqCommands
from scaffolding import PBJelly
import logging
from utils import DirUtils

class RunAssembly(object):
    """
    The RunAssembly class is the main flow of the full pipeline. This class makes the decisions which programs 
    will be executed.
    """
    
    def __init__(self):
        """
        The constructor of the RunAssembly class contains the main flow of the pipeline.
        """
        self.configureLogging()
        logging.info("Starting pipeline")
        pool = self.parseArguments()
        for lib in pool.libs:
            self.preprocess(lib)        
                    
        self.doGenomeSizeEstimation(pool)
        
        self.doAssembly(pool)   
        self.doScaffolding(pool)
        
#         self.contaminationCheck(self.assembly, pool.outputDir)
        self.createStatistics(pool)
        
    def configureLogging(self):
        """
        Setting up the variables for the logging module.
        """
        logging.basicConfig(filename="pipeline.log", format="%(asctime)-25s%(levelname)-12s%(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p")
        #logging.basicConfig(filename="C:/Users/Jetse/workspace/assembly/assembly.log", format="%(asctime)-22s%(levelname)-12s%(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p")
        
        console=logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(logging.Formatter("%(asctime)-25s%(message)s", "%I:%M:%S %p"))
        logging.getLogger().addHandler(console)
            
    def createStatistics(self, pool):
        """
        The method createStatistics creates an assemblystatistics object for generating the statistics of the assembly and creates a fastq report for the raw and fully preprocessed data.
        """
        logging.info("Creating assembly stats")
        assemblyController = AssemblyControl.AssemblyStatistics()
        assemblyController.AssemblyStatisticsOfPipeline(pool.outputDir + "statistics/", pool, self.assembly)
        LaTeX.ltxPart("Supplementary materials")
        for lib in pool.libs:
            self.createFastqReport("raw " + lib.libName, lib.rawForward, lib.rawReversed, lib.outputDir + "raw_qc/")
            self.createFastqReport("preprocessed " + lib.libName, lib.forward, lib.reversed, lib.outputDir + "preprocessed/")
        
        Reporter.instance.createReport(pool.outputDir + "report/")
    
    def doScaffolding(self, pool):
        """
        The method doScaffolding regulates the scaffolding. If pacbio data is available, the scaffolding will be done with pbjelly.
        #TODO add SCARPA scaffolding
        """
        for lib in pool.libs:
            if lib.sequencingPlatform == "pacbio":
                self.assembly = PBJelly.PBJelly(pool.outputDir + "scaffolds/", assembly=self.assembly, reads=lib.forward).execute()
    
    def preprocess(self, lib):
        """
        The method preprocess determines which preprocessing steps have to be executed for a given library.
        """
        logging.info("Preprocessing of " + lib.libName)
        LaTeX.ltxSection("Preprocessing of " + lib.libName)
        lib.forward = DirUtils.fileRegexToList(lib.forward)
        lib.reversed = DirUtils.fileRegexToList(lib.reversed)
        if lib.format == "sff":
            for idx, sffFile in enumerate(lib.forward):
                lib.forward[idx] = SffToFastqConverter.SffToFastqConverter(lib.outputDir, sffFile=sffFile).execute()
        smallReport = FastqSmallReport.FastqSmallReport()
        smallReport.createSmallReport(lib.forward, lib.reversed) 
         
        if len(lib.forward)  > 1:     
            lib.forward = FastqCommands.MergeCommand(lib.outputDir, direction="forward", fastqFiles=lib.forward).execute()
            if lib.reversed != None:
                lib.reversed = FastqCommands.MergeCommand(lib.outputDir, direction="reversed", fastqFiles=lib.forward).execute()
            FastqSmallReport.FastqSmallReport().createSmallReport(lib.forward, lib.reversed)
        else:
            lib.forward = lib.forward[0]
            if lib.reversed != None:
                lib.reversed = lib.reversed[0]
        lib.avgReadlength = float(smallReport.fastqInfo[smallReport.fastqInfo.keys()[0]][2])
        
        if lib.sequencingPlatform == "illumina":
            self.illuminaPreprocess(lib)
        elif lib.sequencingPlatform == "454":
            lib.forward = FastqMcfTrimming.FastqTrimmer(lib.outputDir, forward=lib.forward,noTrim=True).execute()
            FastqSmallReport.FastqSmallReport().createSmallReport(lib.forward, lib.reversed)
        
        self.filterContamination(lib) 
        
    def filterContamination(self, lib):
        """
        The method filterContamination creates a contamination filter for each organism to filter for.
        """
        if lib.sequencingPlatform == "pacbio":
            return
        contamination = [os.path.dirname(__file__) + "/preprocessing/phixDb/PhiX.fasta"]
        if Configuration.instance.getGlobalOption("contamination") != None:
            contamination += Configuration.instance.getGlobalOption("contamination").split(",")
        for cont in contamination:
            [lib.forward,lib.reversed] = ContaminationFiltering.ContaminationFiltering().filterContamination(lib.outputDir, lib.forward, lib.reversed, lib.insertSize, cont)
            FastqSmallReport.FastqSmallReport().createSmallReport(lib.forward, lib.reversed)
        
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
        LaTeX.ltxSection("Assembly")
        if Configuration.instance.getGlobalOption("assembler") == None or Configuration.instance.getGlobalOption("assembler") == "wgs":
            assembler = WgsAssembler.WgsAssembler()
            self.assembly = assembler.doAssembly(pool.outputDir + "assembly/", pool)
        elif Configuration.instance.getGlobalOption("assembler") == "allpaths":
            self.assembly = AllpathsAssembler.AllpathsAssembler().doAssembly(pool.outputDir + "allPathsAssembly/", pool)
            
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
        
        QuakeCorrection.QuakeCorrection().doCorrection(pool, cutoff=genomeSizeEstimator.cutoff)    
        
        
    def illuminaPreprocess(self, lib):
        """
        This method manages the illumina preprocessing by creating fastqMcfTrimming objects to do the trimming/adapter
        removal/duplicate removal.
        """
        #Run fastqmcf to trim the reads
        logging.info("Trimming " + lib.libName)
        if lib.rawReversed !=None:
            fastqTrimmer = FastqMcfTrimming.FastqTrimmer(lib.outputDir, forward=lib.rawForward,reversed=lib.rawReversed)
        else:
            fastqTrimmer = FastqMcfTrimming.FastqTrimmer(lib.outputDir, forward=lib.rawForward)
        lib.forward = fastqTrimmer.execute()
        if lib.rawReversed !=None:
            lib.reversed = fastqTrimmer.outReversed
        FastqSmallReport.FastqSmallReport().createSmallReport(lib.forward, lib.reversed)
        
    
    def createFastqReport(self, libName, forward, reversedFile, outDir):
        """
        This method creates the fastqc report object which executes fastqc.
        """
        logging.info("Creating fastq report for " + libName)
        fastqController = FastqControl.FastqReportCreator(outDir, forward=forward, status="forward reads of " + libName)
        forwardQcIndex = fastqController.execute()
        if reversedFile !=None:
            rawFastqControllerRev = FastqControl.FastqReportCreator(outDir, forward=reversedFile, status="reversed reads of " + libName)
            reversedQcIndex = rawFastqControllerRev.execute()
            return[forwardQcIndex,reversedQcIndex]
        return[forwardQcIndex]
        
if __name__ == '__main__':
    RunAssembly()