import Command, os

class A50Plotter(Command.Command):
        
    def setCommand(self):
        self.outputFile=self.outputDir + "a50Plot.png"
        self.addArg("Rscript")
        self.addArg(os.path.dirname(__file__) +  "/a50plot.R")
        self.addArg(self.faFile)
        self.addArg(self.outputFile)

