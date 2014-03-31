from Bio import SeqIO
from qualityControl.Reporter import LaTeX
import os, logging

class FastqSmallReport(object):
    """
    The clas FastqSmallReport creates a small report with some base statistics for a fastq file. The statistics are:
    * Total bp
    * Total reads
    * Average read length
    * Percentage high quality reads (qv > 30)
    """
    def createSmallReport(self, fastqFiles, libName):
        self.libName = libName
        self.fastqInfo = {}
        for fastqFile in fastqFiles:
            index = os.path.basename(os.path.splitext(fastqFile)[0])
            self.fastqInfo[index] = self.getFastqInfo(fastqFile)
            logging.info(fastqFile + ":")
            logging.info("Total bp: " + self.fastqInfo[index][0])
            logging.info("Total reads: " + self.fastqInfo[index][1])
            logging.info("Avg read length: " + self.fastqInfo[index][2])
            logging.info("percentage good quality: " + self.fastqInfo[index][3])
            
    def getFastqInfo(self, fastqFile):
        totalLength = 0
        totalReads = 0
        totalHighq = 0
        for record in SeqIO.parse(open(fastqFile), "fastq"):
            seqLength = len(record.seq)
            totalReads += 1
            totalLength += seqLength
            for qual in record.letter_annotations["phred_quality"]:
                if qual > 30:
                    totalHighq +=1
            
        return ["{:,}".format(totalLength), "{:,}".format(totalReads), str(totalLength/totalReads),str(int(round(totalHighq/float(totalLength)*100))) + "\%"]
        
    def getLaTeXReport(self):
        txt = ""
        if self.libName != None:
            txt = txt + "\\subsection{Statistics of "+self.libName+"}\n"
        table = LaTeX.ltxTable(len(self.fastqInfo.values()[0])+1)
        table.addRow(["Fastq file","Total bases", "number of reads","AVG read length","Percentage high quality bases ($qv > 30$)"])
        for [index, fastqEntry] in self.fastqInfo.iteritems():
            table.addRow([index] + fastqEntry)
        txt = txt + table.getText()
        return txt
        