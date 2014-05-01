from Bio import SeqIO
from qualityControl.Reporter import LaTeX, Reporter
import os, logging

class FastqSmallReport(object):
    """
    The clas FastqSmallReport creates a small report with some base statistics for a fastq file. The statistics are:
    * Total bp
    * Total reads
    * Average read length
    * Percentage high quality reads (qv > 30)
    """
    def createSmallReport(self, forwardFastq, reversedFastq=None):
        Reporter.instance.objects.append(self)
        self.fastqInfo = {}
        fastqFiles = []
        if type(forwardFastq) is list:
            if reversedFastq == None:
                fastqFiles = forwardFastq
            else:
                fastqFiles = forwardFastq + reversedFastq
        else:
            if reversedFastq == None:
                fastqFiles = [forwardFastq]
            else:
                fastqFiles = [forwardFastq, reversedFastq]
                
        
        for fastqFile in fastqFiles:
            if len(fastqFile) == 0:
                continue
            index = os.path.basename(os.path.splitext(fastqFile)[0])
            self.fastqInfo[index] = self.getFastqInfo(fastqFile)
#             logging.debug(os.path.basename(fastqFile) + ":")
#             logging.debug("Total bp: " + self.fastqInfo[index][0])
#             logging.debug("Total reads: " + self.fastqInfo[index][1])
#             logging.debug("Avg read length: " + self.fastqInfo[index][2])
#             logging.debug("percentage good quality: " + self.fastqInfo[index][3])

            
    def getFastqInfo(self, fastqFile):
        """
        The method getFastqInfo counts the number of reads, bases and high quality reads and writes these to a summary file. 
        If this summary file already exists, this file is used instead of counting the reads/bases again.
        """
        summaryFile = fastqFile + ".summary"
        if os.path.exists(summaryFile):
            return self.getFastqInfoFromSummary(summaryFile)
        logging.info("Counting number of reads of " + fastqFile)
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
        info = ["{:,}".format(totalLength), "{:,}".format(totalReads), str(totalLength/totalReads),str(int(round(totalHighq/float(totalLength)*100))) + "\%"]
        with open(summaryFile, "w") as summaryWriter:
            summaryWriter.write("\t".join(info))  
        return info
    
    def getFastqInfoFromSummary(self, summaryFile):
        """
        This method reads the summary file and returns the content
        """
        with open(summaryFile) as summaryReader:
            for line in summaryReader:
                return line.split("\t")
                
        
    def getLaTeXReport(self):
        txt = ""
        table = LaTeX.ltxTable(len(self.fastqInfo.values()[0])+1)
        table.addRow(["Fastq file","Total bases", "number of reads","AVG read length","Percentage high quality bases ($qv > 30$)"])
        for [index, fastqEntry] in self.fastqInfo.iteritems():
            table.addRow([index] + fastqEntry)
        txt = txt + table.getText()
        return txt + "\\\\"
        