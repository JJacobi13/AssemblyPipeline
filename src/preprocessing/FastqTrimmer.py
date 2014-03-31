from model import Command, Pool
import os, sys

class FastqTrimmer(Command.Command):
    """
    The FastqTrimmer object uses fastq-mcf to trim the adapters from a fastq. Fastq-mcf also removes bad quality ends and duplicate reads.
    The adapters of the fastq file are extracted from the fastqc report.
    """
    
    def checkInput(self):
        if hasattr(self.collection, "forward") == False:
            raise ValueError("Library has no forward reads, which are required for creating a fastqc report")
        
    def setCommand(self):
        self.outputFile = self.collection.outputDir + self.collection.libName + "_trimmed_1.fastq"
        if hasattr(self.collection, "reversed"):
            self.outReversed = self.collection.outputDir + self.collection.libName + "_trimmed_2.fastq"
        
        self.addArg("fastq-mcf")
        self.addArg("-D  60") #TODO: find nice parameter for duplicates, something like avg readlen/2
        self.addArg("-o " + self.outputFile)
        if hasattr(self.collection, "reversed"):
            self.addArg("-o " + self.outReversed)
        adapterFile = self.findAdaptersInQcReport(self.collection)
        self.addArg(adapterFile)
        self.addArg(self.collection.forward)
        if hasattr(self.collection, "reversed"):
            self.addArg(self.collection.reversed)
        
    def updateStatus(self):
        self.collection.status="preprocessed"
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


if __name__ == '__main__':
    if len(sys.argv)==3:
        pool = Pool.Pool(outputDir=sys.argv[2] + "/")
        pool.forward = sys.argv[1]
        pool.fastqcDir = os.path.dirname(sys.argv[1])   
        pool.libName =  os.path.splitext(os.path.basename(sys.argv[1]))[0]
        FastqTrimmer(pool).execute()
    else:
        print("Usage: python FastqControl.py <input fastq> <output directory>")