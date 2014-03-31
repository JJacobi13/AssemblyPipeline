import os
from commandLineCommands import Command

class LatexCreator(Command.Command):
    
    def setCommand(self):
        self.outputFile=os.path.splitext(self.texFile)[0] + ".xpdf"
        self.addArg("pdflatex -interaction=nonstopmode")
        self.addArg("-output-directory="+self.outputDir)
        self.addArg(self.texFile)
