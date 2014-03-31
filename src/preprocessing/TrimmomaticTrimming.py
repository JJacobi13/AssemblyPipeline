from qualityControl import FastqControl
import os
from commandLineCommands import Command
#/mnt/nexenta/haars001/projects/QA_pipeline/filter_hiseq_pipeline.sh
class Trimmomatic(Command.Command):
    """
    Do preprocessing with trimmomatic.
    ERROR: Not updated for pipeline v0.1!
    """    
    def checkInput(self):
        if hasattr(self.collection, "forward") == False:
            raise ValueError("Library has no forward reads, which are required for creating a fastqc report")
        if hasattr(self.collection, "fastqcDir") == False:
            FastqControl.FastqReportCreator(self.collection).execute()
    
    def setCommand(self):
        self.outputFile = self.collection.outputDir + self.collection.libName + "_trimmed_1.fastq"
        if hasattr(self.collection, "reversed"):
            self.outReversed = self.collection.outputDir + self.collection.libName + "_trimmed_2.fastq"
        
        self.addArg("java -jar")
        self.addArg(self.collection.config.getPath("trimmomatic"))
        if hasattr(self.collection, "reversed"):
            self.addArg("PE")   
        else:
            self.addArg("SE")
        self.addArg("-threads " + self.collection.config.getGlobalOption("maxThreads"))
        self.addArg(self.collection.forward)
        if hasattr(self.collection, "reversed"):
            self.addArg(self.collection.reversed)
        self.addArg(self.outputFile)
        self.addArg(self.collection.outputDir + self.collection.libName + "_trimmed_unpaired_1.fastq")
        if hasattr(self.collection, "reversed"):
            self.addArg(self.outReversed)
            self.addArg(self.collection.outputDir + self.collection.libName + "_trimmed_unpaired_2.fastq")
        self.addArg("ILLUMINACLIP:" + self.findAdaptersInQcReport(self.collection) + "2:30:10")
        self.addArg("LEADING:3")
        self.addArg("SLIDINGWINDOW:4:15")
        self.addArg("MINLEN:36")      
        
    def updateStatus(self):
        self.collection.status="trimmomatic"
        self.collection.forward = self.outputFile
        if hasattr(self.collection, "reversed"):
            self.collection.reversed = self.outReversed
    
    def findAdaptersInQcReport(self, lib):
        """
        The method findAdaptersInQcReport extracts the adapters from a fastqc report and writes these to a fasta file.
        #TODO: extract adapters from file if adapters are known
        :param lib: the library to find the adapters for
        :fileType lib: :py:class:`Library.Library`
        :return: a filename of the fasta file all adapters are written to
        """
        forwardReport = lib.fastqcDir + "/" + os.path.splitext(os.path.basename(lib.forward))[0] + "_fastqc/fastqc_data.txt"
        if hasattr(self.collection, "reversed"):
            reversedReport = lib.fastqcDir + os.path.splitext(os.path.basename(lib.reversed))[0] + "_fastqc/fastqc_data.txt"
        adapters = self.findAdapters(forwardReport)
        if hasattr(self.collection, "reversed"):
            adapters += self.findAdapters(reversedReport)
        adapterFile = lib.outputDir+"adapters.fa"
        with open(adapterFile, "w") as fastaWriter:
            for i, adapter in enumerate(adapters):
                fastaWriter.write(">adapter_" + str(i) + "\n")
                fastaWriter.write(adapter + "\n")
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
