from configuration import Configuration
from commandLineCommands import Command
import sys, os, logging, re
from qualityControl.Reporter import LaTeX
from xml.dom import minidom
# from qualityControl import FastqSmallReport

class FastqReportCreator(Command.Command):
    
    def checkInput(self):
        if hasattr(self, "forward") == False:
            raise ValueError("FastqReportCreator has no forward reads, which are required for creating a fastqc report")
    
    def getLaTeXReport(self):
        txt = "\\section{Fastqc of " +self.status+"}\n"
        table = False
        imgs = 0
        with open(self.outputFile) as reportReader:
            for line in reportReader:
                if "</h2>" in line and "Summary" not in line:
                    txt = txt + "\\subsection*{"+re.findall(r"]\">(.*?)</h2>",line)[0]+"}\n"
#                     if "Basic Statistics" in line:
#                         smallReport = FastqSmallReport.FastqSmallReport()
#                         smallReport.createSmallReport([self.forward], None)
#                         txt = txt + smallReport.getLaTeXReport()
#                         txt = txt + "\\\\"
                elif "<table>" in line:
                    xml = "<table>"
                    table = True
                elif "</table>" in line:
                    xml = xml + "</table>"
                    table = False
                    domTable = minidom.parseString(xml)
                    ltxTable = LaTeX.ltxTable(len(domTable.firstChild.firstChild.childNodes))
                    noOfRows = 0
                    for row in domTable.firstChild.childNodes:
                        cols = []
                        noOfRows = noOfRows + 1
                        if noOfRows > 15:
                            continue
                        for col in row.childNodes:
                            cols.append(col.firstChild.nodeValue.replace("%"," percent"))
                        if "Filename" in cols[0]:
                                continue
                        ltxTable.addRow(cols)
                    txt = txt + ltxTable.getText()
                    if noOfRows > 15:
                        txt = txt + "\\\\Total length of this table is "+ str(noOfRows) + ". The table is cut after 15 rows..."
                elif table == True:
                    xml = xml + line.strip()
                elif "<img class=\"indented\"" in line:
                    imgs = imgs +1
                    img = re.findall("src=\"(.*)\" alt=",line)[0]
                    ltxImg = LaTeX.ltxImage(os.path.dirname(self.outputFile) + "/" + img)
                    txt = txt + ltxImg.getText()
                    if imgs % 2 == 0:
                        txt = txt + "\\clearpage\n"
        return txt
                

       
#     def getHtmlReport(self, reportDir):
#         html = "<div class = \"fastqReport\">"
#         html = html + "<h1 id = \""+self.collection.libName+self.libStatus+"\">Fastq report of " + self.libStatus + " reads of library " + self.collection.libName + "</h1>"
#         html = html + "<h2>Forward reads</h2>"
#         fastqReportSummary = self.extractDiv("module", self.forwardReport, reportDir)
#         html = html + fastqReportSummary
#         if hasattr(self, "reversedReport"):
#             html = html + "<h2>Reversed reads</h2>"
#             fastqReportSummary = self.extractDiv("module", self.reversedReport, reportDir)
#             html = html + fastqReportSummary
#          
#         html = html + "</div>"
#         return html
#     
#     def extractDiv(self, divName, report,reportDir):
#         appendToString = 0;
#         divString = ""
#         with(open(report)) as fastqReportReader:
#             for line in fastqReportReader:
#                 if line.startswith("<div class=\""+divName+"\">"):
#                     appendToString = 1
#                 elif line.startswith("</div>") and appendToString == 1:
#                     divString = divString + "<br /><a href=\"" + os.path.relpath(report, reportDir)+"\">full report</a>"
#                     return divString
#                 elif appendToString == 1:
#                     divString = divString + line
    
    
    def setCommand(self):
        self.outDir = self.outputDir
        self.outputFile = self.outDir+os.path.splitext(os.path.basename(self.forward))[0] + "_fastqc/fastqc_report.html"
        
        if not os.path.isdir(self.outDir):
            os.makedirs(self.outDir)
            
        self.addArg("fastqc")
        self.addArg("-o " + self.outDir)
        self.addArg("-t " + Configuration.instance.getGlobalOption("maxThreads"))
        self.addArg(self.forward)
        if hasattr(self, "reversed"):
            self.addArg(self.reversed)
    
if __name__ == '__main__':
    Configuration.instance.setOption("maxThreads", "20")
    Configuration.instance.setOption("overwrite", "0")
    logging.basicConfig(format="%(asctime)-25s%(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p")
    if len(sys.argv) == 3:
        FastqReportCreator(sys.argv[1],forward=sys.argv[2]).execute()
    elif len(sys.argv)==4:
        FastqReportCreator(sys.argv[1],forward=sys.argv[2], reversed=sys.argv[3]).execute()
    else:
        print("Usage: python FastqControl.py <output directory> <forward fastq> [reversed fastq]")
