from model import Command
from Bio import SeqIO
import os, shutil, threading

class Maker(object):
    """
    
    """

    def __init__(self, pool):
        self.pool = pool        
        self.makerDir = self.pool.outputDir + "maker/"
        if not os.path.isdir(self.makerDir):
            os.makedirs(self.makerDir)
        self.changeDefaultMakerOpts()
        
    
    def execute(self):
        fastaFiles = self.splitFasta()
        print len(fastaFiles)
        exit(1)
        threads = []
        for faFile in fastaFiles:
            subDir = self.prepareSubMakerDir(faFile)
            makerCommand = MakerCommand(self.pool, subDir, faFile, self.makerDir)
            threads.append(makerCommand)
            makerCommand.start()
        for t in threads:
            t.join()
        self.mergeOutput(fastaFiles)
        
    def getHtmlReport(self, outDir):
        html = "<div class=\"maker_annotation\">"
        html = html + "<h1>Maker annotation</h1>"
        html = html + "<table>"
        html = html + "<tr><td>Total proteins:</td><td>"+str(self.getNoOfProteins())+"</td></tr>"
        html = html + "</table>"
        html = html + "<h2>Used data</h2>"
        html = html + "<table>"
        html = html + "<tr><td>Est source:</td><td>"+self.pool.config.getMakerOpt("est_source")+"</td></tr>"
        html = html + "<tr><td>Protein source:</td><td>"+self.pool.config.getMakerOpt("protein_source")+"</td></tr>"
        html = html + "<tr><td>Snap source:</td><td>"+self.pool.config.getMakerOpt("snap_source")+"</td></tr>"
        html = html + "<tr><td>Augustus species:</td><td>"+self.pool.config.getMakerOpt("augustus_species").replace(",","<br />")+"</td></tr>"
        html = html + "</table>"
        html = html + "</div>"
        return html
    
    def getNoOfProteins(self):
        records = list(SeqIO.parse(self.makerDir+"allProteins.fasta", "fasta"))
        return len(records)
    
    def mergeOutput(self, fastaFiles):
        gffFiles = []
        faFiles = []
        for faFile in fastaFiles:
            subDir = os.path.dirname(faFile)+"/"+os.path.splitext(os.path.basename(faFile))[0]
            dataStore = subDir + "annotation_master_datastore_index.log"
            
            gffFile = subDir + "/proteins.gff"
            gffFiles.append(gffFile)
            GffMergeCommand(self.pool, outputFile=gffFile, dataStore=dataStore)
            
            faFile = subDir + "/proteins.fasta"
            faFiles.append(faFile)
            FastaMergeCommand(self.pool, outputFile=gffFile, dataStore=dataStore)
            
        GffMergeCommand(self.pool, outputFile=self.makerDir+"allProteins.gff", gffFiles=gffFiles)
        FastaMergeCommand(self.pool, outputFile=self.makerDir+"allProteins.fasta", faFiles=faFiles)
     
    def prepareSubMakerDir(self, faFile):
        subDir = os.path.dirname(faFile)+"/"+os.path.splitext(os.path.basename(faFile))[0]
        if not os.path.isdir(subDir):
            os.makedirs(subDir)
        with open(self.makerDir+"maker_opts.ctl") as optsReader:
            with open(subDir+"/maker_opts.ctl", "w") as optsWriter:
                for line in optsReader:
                    if line.startswith("genome"):
                        line = line.replace("$genome", faFile)
                    optsWriter.write(line)
        return subDir
        
    def splitFasta(self):
        origFasta = self.pool.contigsFasta
        recordIter = SeqIO.parse(open(origFasta),"fasta")
        subFastas = []
        i = 0
        fileNo = 0
        subFastas.append(self.makerDir+"subseq_" + str(fileNo) + ".fasta")
        with open(subFastas[i], "w") as fastaWriter:
            for read in recordIter:
                if i > fileNo:
                    fastaWriter.close()
                    fileNo = fileNo +1
                    subFastas.append(self.makerDir+"subseq_" + str(fileNo) + ".fasta")
                    fastaWriter = open(subFastas[i], "w")
                    i = 0
                i = i+1
                    
                SeqIO.write(read, fastaWriter, "fasta")
        return subFastas       
        
    def changeDefaultMakerOpts(self):
        scriptDir = os.path.dirname(os.path.realpath(__file__)) + "/"
        optsFile = self.makerDir+"maker_opts.ctl"
        shutil.copyfile(scriptDir+"maker_bopts.ctl", self.makerDir+"maker_bopts.ctl")
        shutil.copyfile(scriptDir+"maker_exe.ctl", self.makerDir+"maker_exe.ctl")
        with open(scriptDir+"maker_opts.ctl") as optsReader:
            with open(optsFile, "w") as optsWriter:
                for line in optsReader:
                    if line.startswith("organism_type"):
                        line = line.replace("$orgType", self.pool.config.getMakerOpt("organism_type"))
                    elif line.startswith("est"):
                        line = line.replace("$est", self.pool.config.getMakerOpt("est"))
                    elif line.startswith("protein"):
                        line = line.replace("$protein", self.pool.config.getMakerOpt("protein"))
                    elif line.startswith("snaphmm"):
                        line = line.replace("$snapHmm", self.pool.config.getMakerOpt("snaphmm"))
                    elif line.startswith("augustus_species"):
                        line = line.replace("$augustusSpecies", self.pool.config.getMakerOpt("augustus_species"))
                    optsWriter.write(line)
                    
class FastaMergeCommand(Command.Command):
    def __init__(self, pool,outputFile,  faFiles=None, dataStore=None):
        self.faFiles = faFiles
        self.outputFile = outputFile
        self.dataStore = dataStore
        Command.Command.__init__(self,pool)
        
    def setCommand(self):
        self.addArg("fasta_merge")
        if self.faFiles !=None:
            self.addArg("-i")
            for faFile in self.faFiles:
                self.addArg(faFile)
        elif self.dataStore != None:
            self.addArg("-d " + self.dataStore)
        self.addArg("-o " + self.outputFile)
        
    def updateStatus(self):
        pass
                      
class GffMergeCommand(Command.Command):
    def __init__(self, pool, outputFile, gffFiles=None, dataStore=None):
        self.gffFiles = gffFiles
        self.outputFile = outputFile
        self.dataStore = dataStore
        Command.Command.__init__(self,pool)
        
    def setCommand(self):
        self.addArg("gff3_merge")
        if self.gffFiles !=None:
            for gffFile in self.gffFiles:
                self.addArg(gffFile)
        elif self.dataStore != None:
            self.addArg("-d " + self.dataStore)
        self.addArg("-o " + self.outputFile)
        
    def updateStatus(self):
        pass
        
class MakerCommand(Command.Command, threading.Thread):
    def __init__(self, pool, subDir, fastaFile, outDir):
        self.subDir = subDir
        self.fastaFile = fastaFile
        self.outDir = outDir
        threading.Thread.__init__(self)
        Command.Command.__init__(self,pool)
    
    def run(self):
        self.execute()
        
    def setCommand(self):
        self.addArg("maker")
        self.addArg("-base " + self.subDir + "/annotation")
        self.addArg(self.subDir+"/maker_opts.ctl")
        self.addArg(self.outDir+"maker_bopts.ctl")
        self.addArg(self.outDir+"maker_exe.ctl")
    
    def updateStatus(self):
        pass