import logging, os, subprocess, re, argparse
from Bio import SeqIO
from commandLineCommands import Rscripts, Mappers, BamCommands, FastaCommands
from qualityControl.Reporter import Reporter, LaTeX
from configuration import Configuration
# from qualityControl import BlastScanner

class AssemblyStatistics(object):
    """
    The class AssemblyStatistics calculates the following statistics of a given fasta file:
    * Total sequences
    * Total length
    * GC percentage
    * Longest sequence
    * N50
    * N90
    * Cegma (percentage of the core eukaryotic genes)
    * Heterozygosity
    * Dna mapping percentage
    * Rna mapping percentage
    * Generates an a50 plot
    
    All statistics and errors are written to the log and a pfd file is created.
    
    When this class is executed from the pipeline, all variables are extracted from the pool and library objects, else those are read from the commandline. 
    These variables are used for generating the statistics.
    """
    
    def main(self):
        """
        When only statistics of an fasta file have to be calculated, call this method. 
        This method parses the commandline arguments and executes the quality control of the assembly.
        """
        logging.basicConfig(format="%(asctime)-25s%(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p")
        
        parser = argparse.ArgumentParser(description="Create statistics of a generated assembly")
        parser.add_argument("fastaFile", type=argparse.FileType('r'), help="The generated assembly in fasta format")
        parser.add_argument("outDir", help="The ouput directory")
        parser.add_argument("--forwardReads", help="The raw forward reads", required=False)
        parser.add_argument("--reversedReads", help="The raw reversed reads", required=False)
        parser.add_argument("--insertSize", help="The insert size between the forward and reversed reads", required=False, metavar="N")
        parser.add_argument("--forwardRna", help="The raw forward RNA reads", required=False)
        parser.add_argument("--reversedRna", help="The raw reversed RNA reads", required=False)
        parser.add_argument("--nrDb", help="an already downloaded nr database for contamination check", required=False)
        parser.add_argument("--maxThreads", help="The number of threads to use (default 20)", metavar="N", default=20)
        parser.add_argument("--overwrite", help="Overwrite existing files", action="store_true")
        args = parser.parse_args()
        args.fastaFile.close()
        self.fastaFile = args.fastaFile.name
        
        
        self.outputDir = args.outDir + "/"
        self.forward = args.forwardReads
        self.reversed = args.reversedReads
        self.insertSize = args.insertSize
        self.forwardRna = args.forwardRna
        self.reversedRna = args.reversedRna
        self.nrDb = args.nrDb
        
        Configuration.instance.setOption("maxThreads", str(args.maxThreads))
        if args.overwrite == True:
            Configuration.instance.setOption("overwrite", "1")
        else:
            Configuration.instance.setOption("overwrite", "0")
        
        self.doQualityControl()
        Reporter.instance.createReport(self.outputDir)
    
    def AssemblyStatisticsOfPipeline(self, outputDir, pool, fastaFile):
        """
        When the statistics of the assembly of the pipeline have to be calculated, call this method. 
        This method parses the Pool and library objects to instance variables and executes the quality control of the assembly.
        """
        self.fastaFile = fastaFile
        self.outputDir = outputDir
        self.forward = pool.libs[0].rawForward
        self.reversed = pool.libs[0].rawReversed
        self.insertSize = pool.libs[0].insertSize
        if hasattr(pool, "forwardRna"):
            self.forwardRna = pool.forwardRna
        else:
            self.forwardRna = None
            
        if hasattr(pool, "reversedRna"):
            self.reversedRna = pool.reversedRna
        else:
            self.reversedRna = None
        self.nrDb = Configuration.instance.getGlobalOption("nrDb")
        self.doQualityControl()
    
    def getLaTeXReport(self):
        """
        Convert all previously calculated statistics into LaTeX with this method.
        """
        txt = "\\section{Assembly statistics}\n"
        table = LaTeX.ltxTable(2)
        table.addRow(["Total sequences: ",str(self.totalSeqs)])
        table.addRow(["Total length: ","{:,}".format(self.totalLen)])
        table.addRow(["GC perc: ","{:.2f}".format(self.gcPerc) + "\%"])
        table.addRow(["Longest sequence: ","{:,}".format(self.longestSeq)])
        table.addRow(["N50 index: ","{:,}".format(self.n50Index)])
        table.addRow(["N50: ","{:,}".format(self.n50)])
        table.addRow(["",""])
        table.addRow(["N90 index: ","{:,}".format(self.n90Index)])
        table.addRow(["N90: ","{:,}".format(self.n90)])
        if hasattr(self, "cegmaScore"):
            table.addRow(["",""])
            table.addRow(["Cegma complete: ",self.cegmaScore[0] + "\%"])
            table.addRow(["Cegma partial: ",self.cegmaScore[1] + "\%"])
        if hasattr(self, "rawDnaMappingStats"):
            table.addRow(["",""])
            table.addRow(["DNA reads: ","{:,}".format(int(self.rawDnaMappingStats["total"]))])
            table.addRow(["Mapped: ",self.rawDnaMappingStats["mapped"] + "\%"])
            if "propPair" in self.rawDnaMappingStats:
                table.addRow(["Properly paired",self.rawDnaMappingStats["propPair"] + "\%"])
            table.addRow(["SNP density: ","{:.2f}".format(self.snpDensity) + " SNPs per 10kb"])
        if hasattr(self, "rnaMappingStats"):
            table.addRow(["",""])
            table.addRow(["RNA reads: ","{:,}".format(int(self.rnaMappingStats["total"]))])
            table.addRow(["Mapped: ",self.rnaMappingStats["mapped"] + "\%"])
            if "propPair" in self.rnaMappingStats:
                table.addRow(["Properly paired: ",self.rnaMappingStats["propPair"] + "\%"])
            
        txt = txt + table.getText()
        txt = txt + "\\begin{figure}[h]\n"
        txt = txt + "\\includegraphics[scale=0.7]{" + self.a50Plot + "}\n"
        txt = txt + "\\end{figure}\n"
        return txt
       
    def doQualityControl(self): 
        """
        This method is the mainflow of the program, which consists of the following steps:
        * calculate the base statistics
        * calculate the SNP density
        * calculate the RNA mapping percentage
        * create an A 50 plot
        * calculate the cegma score
        """
        Reporter.instance.objects.append(self)
        
        self.calculateBaseStatistics()
        if self.forward != None:
            self.calculateSnpDensity()
        if self.forwardRna != None:
            self.getRnaMappingPerc()
        self.a50Plot = Rscripts.A50Plotter(self.outputDir, faFile = self.fastaFile).execute()
#         self.cegmaScore = self.getCegmaStatistics()      
        
    def getCegmaStatistics(self):
        """
        This method executes the cegma command and parses the output file to retrieve the percentage of core eukaryotic genes which are found complete or partial.
        """
        cegmaFile = FastaCommands.CegmaCommand(self.outputDir, genome=self.fastaFile).execute()
        complete = "-"
        partial = "-"
        with open(cegmaFile) as fileReader:
            for line in fileReader:
                if not line.strip() or line.startswith("#"):
                    continue
                info = line.split()
                if info[0] == "Complete":
                    complete = info[2]
                elif info[0] == "Partial":
                    partial = info[2]
        return [complete,partial]
        
    def getRnaMappingPerc(self):
        """
        This method executes Tophat to map the RNA reads against the assembly, merges the output bam fils of tophat (mapped/unmapped) and calculates the percentage of mapped reads.
        """
        mappedBam = Mappers.Tophat(self.outputDir, refGenome=self.fastaFile, forwardRna=self.forwardRna, reversedRna=self.reversedRna).execute()
        mergedBam = BamCommands.BamMerger(self.outputDir, bamFiles = [mappedBam,os.path.dirname(mappedBam) +"/unmapped.bam"]).execute()
        self.rnaMappingStats = self.getMappingPerc(mergedBam)
        
       
    def calculateSnpDensity(self):
        """
        This method maps the reads with bowtie against the reference genome and converts the output to bam format.
        After this the mapping percentage is calculated and a SNP calling is executed.
        The SnpDensity is defined as the chance a base is a SNP
        """
        if self.reversed == None:
            samFile = Mappers.Bowtie(self.outputDir, refGenome=self.fastaFile, fastqFile=self.forward).execute()
        else:
            samFile = Mappers.Bowtie(self.outputDir, refGenome=self.fastaFile, forward=self.forward, reversed=self.reversed, insertSize=self.insertSize).execute()

        bamFile = BamCommands.SamToBamConverter(self.outputDir, samFile=samFile).execute()
        self.rawDnaMappingStats = self.getMappingPerc(bamFile)
        vcfFile = BamCommands.SamtoolsMpileup(self.outputDir, bamFile=bamFile, fastaFile=self.fastaFile).execute()
        self.snpDensity = self.getSnpDensity(vcfFile)
        logging.info("SNP density: " + "{:.2E}".format(self.snpDensity))
           
    def calculateBaseStatistics(self):
        """
        This method calculates the base statistics which consist of the following statistics:
        * total number of sequences
        * total length
        * GC percentage
        * longest sequence
        * N50
        * N90
        all output is written to the log
        """
        [contigLengths, totalLength] = self.getContigsLengths(self.fastaFile)
        [n50Index, n50] = self.calculateN(50, contigLengths, totalLength)
        [n90Index, n90] = self.calculateN(90, contigLengths, totalLength)
        
        self.totalSeqs = len(contigLengths)
        self.totalLen = totalLength
        self.gcPerc = self.gcNo/float(totalLength)*100
        self.longestSeq = contigLengths[0]
        self.n50 = n50
        self.n50Index = n50Index
        self.n90 = n90
        self.n90Index = n90Index
        
        logging.info("Total sequences: " + str(self.totalSeqs))
        logging.info("Total length: " + str(totalLength))
        logging.info("GC perc: " + str(self.gcPerc) + "%")
        logging.info("Longest sequence: " + str(self.longestSeq))
        logging.info("N50 index: " + str(n50Index) + ", N50: " + str(n50))
        logging.info("N90 index: " + str(n90Index) + ", N90: " + str(n90))
        
    def getSnpDensity(self, vcfFile):
        """
        The SNP density is calculated by reading all lines of a vcf files which are not comments and are not empty
        """
        noOfSnps = 0
        with open(vcfFile) as vcfReader:
            for line in vcfReader:
                if line.startswith("#") or line.strip() == False:
                    continue
                noOfSnps = noOfSnps +1
        return noOfSnps/float(self.totalLen)*10000
    
    def getMappingPerc(self, bamFile):
        """
        The mapping percentage is calculated by samtools flagstat, the output of samtools flagstat is parsed and put into a dictionary
        """
        results = {}
        p = subprocess.Popen("samtools flagstat " + bamFile, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in iter(p.stdout.readline, ''):
            if len(line) > 0:
                logging.debug(line.strip())
            if("mapped (" in line):
                results["mapped"] = self.RnaMappedPerc = re.findall("\d+.\d+", line)[1]
            elif "properly paired" in line and self.reversed != None:
                results["propPair"] = re.findall("\d+.\d+", line)[1]
            elif "in total" in line:
                results["total"] = re.findall("\d+.\d+", line)[0]
        return results
        
    def getContigsLengths(self, fastaFile):
        """
        The lengths of the contigs are calculated with the biopython library. In the end all contigs are sorted on length.
        """
        contigLengths = []
        totalLength = 0
        self.gcNo = 0
        for contig in SeqIO.parse(open(fastaFile), "fasta"):
            self.gcNo = self.gcNo + sum(map(contig.seq.count, ['G', 'C', 'g', 'c']))
            contigLenght = len(contig.seq)
            totalLength += contigLenght
            contigLengths.append(contigLenght)
        
        contigLengths.sort()
        contigLengths.reverse()
        return [contigLengths, totalLength]
    
    def calculateN(self, n, contigLengths, totalLength):
        """
        The method calculateN calculates the N50, N90 or any other N.
        The N50 is defined as
        "Given a set of sequences of varying lengths, the N50 length is defined as the length N for which 50% of all bases in the sequences are in a sequence of length L < N."
        """
        lengthTillNow = 0
        indexTillNow = 0 
        toThisLength = totalLength * n / 100
        for contigLen in contigLengths:
            lengthTillNow += contigLen
            indexTillNow = indexTillNow + 1
            if toThisLength < lengthTillNow:
                return [indexTillNow,contigLen]
               
if __name__ == '__main__':
    AssemblyStatistics().main()