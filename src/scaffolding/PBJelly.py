from commandLineCommands import Command
import os
from configuration import Configuration

class PBJelly(Command.Command):
    """
    The PBJelly object executes all 5 steps of PBJelly in a row.
    #TODO: !!!NOT hardcoded path to Jelly.py!!!
    """    
    def checkInput(self):
        self.createProtocol()
        
    def createProtocol(self):
        self.protocol = self.outputDir + "protocol.xml"
        with open(self.protocol, "w") as protocolWriter:
            protocolWriter.write("<jellyProtocol>")
            protocolWriter.write("<blasr>-minMatch 8 -minPctIdentity 70 -bestn 5 -nCandidates 20 -maxScore -500 -nproc 4 -noSplitSubreads</blasr>")
            protocolWriter.write("<reference>"+self.assembly+"</reference>")
            protocolWriter.write("<outputDir>"+self.outputDir+"</outputDir>")
            protocolWriter.write("<input baseDir=\""+os.path.dirname(self.reads)+"\">")
            protocolWriter.write("<job>"+os.path.basename(self.reads)+"</job>")
            protocolWriter.write("</input>")
            protocolWriter.write("</jellyProtocol>")
    
    def setCommand(self):
        self.outputFile = self.outputDir + "jelly.out.fasta"
        self.addArg("python /home/jaco001/fungiAssembly/programs/Jelly_13.10.22/pbsuite/jelly/Jelly.py setup " + self.protocol + ";")
        self.addArg("python /home/jaco001/fungiAssembly/programs/Jelly_13.10.22/pbsuite/jelly/Jelly.py mapping " + self.protocol + ";")
        self.addArg("python /home/jaco001/fungiAssembly/programs/Jelly_13.10.22/pbsuite/jelly/Jelly.py support " + self.protocol + ";")
        self.addArg("python /home/jaco001/fungiAssembly/programs/Jelly_13.10.22/pbsuite/jelly/Jelly.py extraction " + self.protocol + ";")
        self.addArg("python /home/jaco001/fungiAssembly/programs/Jelly_13.10.22/pbsuite/jelly/Jelly.py assembly " + self.protocol + " -x \"--nproc="+Configuration.instance.getGlobalOption("maxThreads")+"\";")
        self.addArg("python /home/jaco001/fungiAssembly/programs/Jelly_13.10.22/pbsuite/jelly/Jelly.py output " + self.protocol + ";")