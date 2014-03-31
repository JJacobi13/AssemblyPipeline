from commandLineCommands import Command
from model import Pool
from qualityControl import FastqControl
import os, sys

class FastqTrimmer(Command.Command):
    """
    The FastqTrimmer object uses fastq-mcf to trim the adapters from a fastq. Fastq-mcf also removes bad quality ends and duplicate reads.
    The adapters of the fastq file are extracted from the fastqc report.
    Parameters:
    * forward
    * reversed
    * fastqcDir
    """               
    def setCommand(self):
        self.outputFile = self.outputDir + "fastqMcfTrimmed_1.fastq"
        if hasattr(self, "reversed"):
            self.outReversed = self.outputDir + "fastqMcfTrimmed_2.fastq"
        
        self.addArg("fastq-mcf")
        self.addArg("-D  60") #TODO: find nice parameter for duplicates, something like avg readlen/2
        if hasattr(self, "noTrim"):
            self.addArg("-0")
        self.addArg("-o " + self.outputFile)
        if hasattr(self, "reversed"):
            self.addArg("-o " + self.outReversed)
        adapterFile = self.findAdaptersInQcReport()
        self.addArg(adapterFile)
        self.addArg(self.forward)
        if hasattr(self, "reversed"):
            self.addArg(self.reversed)
    
    def findAdaptersInQcReport(self):
        """
        The method findAdaptersInQcReport extracts the adapters from a fastqc report and writes these to a fasta file.
        :param lib: the library to find the adapters for
        :fileType lib: :py:class:`Library.Library`
        :return: a filename of the fasta file all adapters are written to
        #TODO: extract adapters from file if adapters are known
        """
        adapters = []
        for fastqReport in self.fastqReports:
            adapters = adapters + self.findAdapters(os.path.dirname(fastqReport) + "/fastqc_data.txt")
        adapterFile = self.outputDir+"adapters.fa"
        with open(adapterFile, "w") as fastaWriter:
            for i, adapter in enumerate(adapters):
                fastaWriter.write(">adapter_" + str(i) + "\n")
                fastaWriter.write(adapter[i] + "\n")
        return adapterFile
        
    def findAdapters(self, report):
        """
        The method findAdapters extracts the adapters of a given fastqc report.
        :return: An array of adapters
        """
        addAdapters=False
        adapters = []
        with open(report) as reportReader:
            for line in reportReader:
                if line.startswith("#"):
                    continue
                elif addAdapters == True:
                    if line.startswith(">>END_MODULE"):
                        return adapters
                    adapters.append(line.split("\t")[0])
                elif line.startswith(">>Overrepresented sequences"):
                    addAdapters = True
