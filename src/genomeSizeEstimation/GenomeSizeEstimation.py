from commandLineCommands import JellyFish
from Bio import SeqIO

import os, logging, subprocess, argparse, locale, sys
from configuration import Configuration
from model import Pool,Library
from argparse import _ensure_value
from qualityControl.Reporter import LaTeX, Reporter

class GenomeSizeEstimation(object):    
    """
    This object is based on martijn Derks is genome estimation scripts    
    """
    def doGenomeSizeEstimation(self, outputDir, pool):
        """
        The method doGenomeSizeEstimation contains the mainflow of the genomesize estimation. This mainflow contains the following methods:
        * Execute Jellyfish count
        * Execute Jellyfish stats
        * Create a histogram of the unique kmers with Jellyfish histo
        * Draw a histogram of the unique kmers
        * Estimate the genome size with the BGI method
        """
        logging.info("Starting genome size estimation")
        if not os.path.isdir(outputDir):
            os.makedirs(outputDir)
            
        LaTeX.ltxSection("Genome size estimation}")
        self.jellyFishCountsFile = JellyFish.JellyFishCount(outputDir, pool=pool).execute()
        self.jellyFishStatsFile = JellyFish.JellyFishStats(outputDir, jellyFishCountsFile=self.jellyFishCountsFile).execute()
        self.jellyfishHistoFile = JellyFish.JellyFishHisto(outputDir, jellyFishCountsFile=self.jellyFishCountsFile).execute()
        self.genSizeHistoPlot = outputDir + "kmer_graph.png"
        self.peak = int(self.drawHisto(self.jellyfishHistoFile, self.genSizeHistoPlot))
        self.calculateGenomeSize(pool, self.jellyFishStatsFile, self.jellyfishHistoFile)
        Reporter.instance.objects.append(self)
    
    def getLaTeXReport(self):
        table = LaTeX.ltxTable(2)
        table.addRow(["Total bp: ","{:,}".format(self.totalBp)])
        table.addRow(["Peak: ",str(self.peak)])
        table.addRow(["Mean base coverage: ", "{:.2f}".format(self.coverage)])
        table.addRow(["",""])
        table.addRow(["Unique good kmers: ","{:,}".format(self.unique_gkmers)])
        table.addRow(["BGI genome size estimation: ", "{:,}".format(int(round(self.bgi)))])
        table.addRow(["GSE kmers/peak: ","{:,}".format(self.kmersPerPeak)])
        
        tex = table.getText()
        
        img = LaTeX.ltxImage(self.genSizeHistoPlot)
        tex = tex + img.getText()
        return tex
    
    def getDescription(self):
        return """The genome size estimation is calculated with three different methods. \\\\
The first method is by counting all unique good kmers (all kmers after the first valley of the graph). \\\\
The second method is the BGI method, this is the most used method to calculate the genome size. 
The BGI method has the following formula: $$BGI = total bp/coverage$$ 
The coverage is calculated with the following formula: 
$$coverage=peak/\\left(\\frac{mean readlength - kmer size + 1}{mean readlength}\\right)$$ 
The peak is calculated by calculating the top of a smooth histogram of the k-mer counts, 
the arrow in the kmer coverage histogram point to the peak. When multiple peaks are found, the first peak is used.\\\\
The third method is dividing the total number of good kmers by the peak.\\\\
For counting the number of kmers and creating a histogram of the kmers jellyfish version 1.1.11 is used (WARNING, hardcoded)\\\\"""
#     
#     def getHtmlReport(self, reportDir):
#         html = "<div class = \"genomesize_estimation\">"
#         html = html + "<h1 id=\"gensize\">Genome size estimation</h1>"
#         html = html + "<table>"
#         html = html + "<tr><td>Total bp:</td><td>"+"{:,}".format(self.totalBp)+"</td></tr>"
#         html = html + "<tr><td>Peak:</td><td>"+str(self.peak)+"</td></tr>"
#         html = html + "<tr><td>Mean base coverage:</td><td>"+"{:.2f}".format(self.coverage)+"</td></tr>"
#         html = html + "<tr><td>BGI genome size estimation:</td><td>"+"{:,}".format(int(round(self.bgi)))+"</td></tr>"
#         html = html + "</table>"
#         html = html + "<img src=\""+os.path.relpath(self.genSizeHistoPlot, reportDir)+"\" />"
#         html = html + "</div>"
#         return html
        
    def drawHisto(self,inputFile, outputFile):
        proc = subprocess.Popen(["Rscript",os.path.dirname(__file__) + "/plot_histo.R",inputFile,outputFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        peak,stErr = proc.communicate()
        logging.debug("draw histo out: " + peak)
        if len(stErr) > 1:
            logging.debug(stErr)
        splittedOut = peak.split()
        if splittedOut[2] != "null":
            logging.warn("Warning: multiple peaks found, using the first peak")
            i = 1
            while splittedOut[i] != "null":
                logging.debug(str(i) + "=" + splittedOut[i])
                i = i+1
        return splittedOut[1] #If there are more peak, return the first, first argument is [1] for first print statement
        
    def calculateGenomeSize(self, pool, jellyFishStatsFile, jellyFishHistoFile):
        """Calculate the following statistics as instance variables:
        * reads: The total number of reads
        * totalBp: The total number of basepairs
        * bkmers: The total of bad kmers (including duplicates)
        * unique_gkmers: The total of unique good kmers
        * totalKmers: The total number of kmers (including duplicates)
        * coverage: 
        """
        self.read_fastq(pool)
        self.pool = pool
        self.calculate_bkmers_gkmers(jellyFishStatsFile, jellyFishHistoFile)
        
        self.coverage = int(self.peak)/((float(self.meanReadLen) - int(Configuration.instance.getGlobalOption("kmer")) + 1)/float(self.meanReadLen))
        self.bgi = self.totalBp/self.coverage
        self.kmersPerPeak = int((int(self.totalKmers)-int(self.bkmers))/float(self.peak))
        
        logging.info("Total bp: "+"{:,}".format(self.totalBp))
        logging.info("Total reads: "+"{:,}".format(self.reads))
        logging.info("meanReadLen: " + str(self.meanReadLen))
        logging.info("Peak: "+str(self.peak))
        logging.info("Mean base coverage: "+ "{:.2f}".format(self.coverage))
        logging.info("Unique good kmers: "+"{:,}".format(self.unique_gkmers))
        logging.info("BGI genome size estimation: "+ "{:,}".format(int(round(self.bgi))))
        logging.info("Total kmers: " + "{:,}".format(self.totalKmers))
        logging.info("Bad kmers: " + "{:,}".format(self.bkmers))
        
        logging.info("GSE kmers/peak: "+"{:,}".format(self.kmersPerPeak))
       
        
    def read_fastq(self, pool):
        self.reads = 0
        self.totalBp = 0
        locale.setlocale( locale.LC_ALL, 'en_US.UTF-8' ) 
        for lib in pool.libs:
            if os.path.exists(lib.forward + ".summary"):
                with open(lib.forward + ".summary") as summaryReader:
                    for line in summaryReader:
                        self.reads += locale.atoi(line.split("\t")[1])
                        self.totalBp += locale.atoi(line.split("\t")[0])
                if lib.reversed !=None:
                    if os.path.exists(lib.reversed + ".summary"):
                        with open(lib.reversed + ".summary") as summaryReader:
                            for line in summaryReader:
                                self.reads += locale.atoi(line.split("\t")[1])
                                self.totalBp += locale.atoi(line.split("\t")[0])
            else:    
                reads = 0
                totalBp = 0
                for rec in SeqIO.parse(lib.forward, "fastq"):
                    reads += 1
                    totalBp += len(rec)
                if lib.reversed != None:
                    self.totalBp = totalBp * 2
                    self.reads = reads * 2
                else:
                    self.totalBp += totalBp
                    self.reads += reads
             
        self.meanReadLen = float(self.totalBp)/self.reads
        
    def get_total_kmers(self, genomeSizeStatsFile): ## Get total kmers from jellyfish stats file
        with open(genomeSizeStatsFile, "r") as stats:
            for line in stats:
                if line.split(":")[0].startswith("Total"):
                    return int(line.split(":")[1])
        
    def calculate_bkmers_gkmers(self, genSizeFile, genSizeHisto): ## Calculate bad and good k-mers
        self.bkmers = 0
        self.cutoff = 0
        self.unique_gkmers = 0
        self.totalKmers = self.get_total_kmers(genSizeFile)
        with open(genSizeHisto, "r") as table:
            prev = self.totalKmers
            bad=True
            for kmer in table:
                kmer = kmer.split(" ")
                if int(kmer[1]) < prev and bad == True:
                    prev = int(kmer[1])
                    self.bkmers += (int(kmer[1]) * int(kmer[0]))
                    self.cutoff = kmer[1]
                else:
                    bad = False
                    self.unique_gkmers += int(kmer[1])
        self.gkmers = int(self.totalKmers) - self.bkmers
 
def required_length_multi(nmin,nmax):
    class RequiredLength(argparse.Action):
         
        def __call__(self, parser, namespace, values, option_string=None):
            if not nmin<=len(values)<=nmax:
                msg='argument "{f}" requires between {nmin} and {nmax} arguments'.format(f=self.dest, nmin=nmin, nmax=nmax)
                raise argparse.ArgumentTypeError(msg)
            items = _ensure_value(namespace, self.dest, [])
            items.append(values)
            setattr(namespace, self.dest, items)
    return RequiredLength
 
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Do a genome size estimation on a set of raw NGS reads')
    parser.add_argument("--outputDir", required=True, help="The output directory to write to")
    parser.add_argument("--lib", required=True, nargs="*", action=required_length_multi(1,2), help="use as: --lib <forward reads> [reversed reads]")
    namespaces = parser.parse_args(sys.argv[1:])

    logging.basicConfig(format="%(asctime)-25s%(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p")
    Configuration.instance.setOption("kmer", "17")
    Configuration.instance.setOption("expCoverage", "50")
    Configuration.instance.setOption("maxMem", "200000000")
    Configuration.instance.setOption("maxThreads", "20")
    Configuration.instance.setOption("overwrite", "0")
     
    #TODO create a pool object, add the libraries and execute the genome size estimation
    pool = Pool.Pool(namespaces.outputDir)
    for lib in namespaces.lib:
        libObj = Library.Library(pool, os.path.splitext(os.path.basename(lib[0]))[0])
        libObj.forward = lib[0]
        if len(lib) > 1:
            libObj.reversed = lib[1]
        pool.libs.append(libObj)
    GenomeSizeEstimation().doGenomeSizeEstimation(pool.outputDir, pool)
