import os
from xml.dom import minidom
from commandLineCommands import Command

class Reporter(object):
    
    def __init__(self, ):
        self.objects = []
        
    def createReport(self, outputDir):
        self.outputDir = outputDir
        if not os.path.isdir(self.outputDir):
            os.makedirs(self.outputDir)
#         self.createHtmlReport(objects)
        self.createLaTeXReport(self.objects)
        
    def createLaTeXReport(self, objects):
        self.fileName = self.outputDir + "report.tex"
        with open(self.fileName, "w") as outWriter:
            outWriter.write("\\documentclass{report}\n")
            outWriter.write("\\usepackage[margin=2cm]{geometry}\n")
            outWriter.write("\\usepackage{graphicx}\n")
#             outWriter.write("\\usepackage{hyperref}\n")
#             outWriter.write("\\hypersetup{colorlinks,citecolor=black, filecolor=black,linkcolor=black,urlcolor=black}\n")
            
            outWriter.write("\\begin{document}\n")
            
            outWriter.write("\\title{Report of the assembly pipeline}\n")
            outWriter.write("\\author{Automatically generated}\n")
            outWriter.write("\\date{\\today}\n")
            outWriter.write("\\maketitle\n")
            
            outWriter.write("\\tableofcontents\n")
            outWriter.write("\\setcounter{secnumdepth}{0}\n")
            
            for objectInstance in objects:
                if hasattr(objectInstance, "getLaTeXReport"):
                    outWriter.write("\\clearpage")
                    outWriter.write(objectInstance.getLaTeXReport())
            outWriter.write("\\end{document}\n")
        
        latexCreator = LatexCreator(self.outputDir, texFile=self.fileName)
        #execute twice, so table of contents is generated too.
        latexCreator.execute()
        latexCreator.execute()

#     def createHtmlReport(self, *objects):
#         self.indexFile = self.outputDir + "index.html"
#         if not os.path.isdir(self.outputDir):
#             os.makedirs(self.outputDir)
#               
#         with open(self.indexFile, "w") as outWriter:
#             outWriter.write("<div class=\"main\">")
#             for objectInstance in objects:
#                 if hasattr(objectInstance,"getHtmlReport") and callable(getattr(objectInstance,"getHtmlReport")):
#                     
#                     outWriter.write(objectInstance.getHtmlReport(self.outputDir))
#             outWriter.write("</div>")
#                     
#         with open(self.indexFile + "new", "w") as outWriter:
#             outWriter.write("<html><head><title>Control of the assembly pipeline</title><link rel=\"stylesheet\" type=\"text/css\" href=\"layout.css"+"\"></head>") 
#             outWriter.write("<div class=\"header\">Control report of the assembly pipeline</div>")
#             shutil.copy(os.path.dirname(__file__) + "/layout.css", self.outputDir + "layout.css")
#             outWriter.write(self.createIndex())
#             with open(self.indexFile) as indexFileReader:
#                 for line in indexFileReader:
#                     outWriter.write(line)
#             outWriter.write("</html>")
#         shutil.move(self.indexFile + "new", self.indexFile)
           
    def createIndex(self):
        html = "<div class=\"summary\">"
        html = html + "<h1>Index</h1>"
        html = html + "<lu>"
        
        doc = minidom.parse(self.indexFile)
        headers = doc.getElementsByTagName("h1")
        for header in headers:
            html = html + "<li><a href=\"#" + header.getAttribute("id") +"\">"+header.firstChild.nodeValue+"</a></li>"
        html = html + "</lu>"
        html = html + "</div>"
        return html

class LatexCreator(Command.Command):
    
    def setCommand(self):
        self.outputFile=os.path.splitext(self.texFile)[0] + ".xpdf"
        self.addArg("pdflatex -interaction=nonstopmode")
        self.addArg("-output-directory="+self.outputDir)
        self.addArg(self.texFile)

instance = Reporter()        
        