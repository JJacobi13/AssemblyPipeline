import os, sys, math, re, logging
from  commandLineCommands import Mappers, Command
from configuration import Configuration
from Reporter import LaTeX, Reporter

class InsertSizeChecker(object):
    """
    The InsertSizeChecker class determines the insert size by mapping with bowtie against a reference genome and compare
    the forward mapping position with the reversed mapping position.
    """
    def main(self):
        """
        The main method is executed when executing the InsertSizeChecker directly from the command line. All arguments are
        retrieved from the command-line.
        """
        if len(sys.argv) != 7:
            print "usage: InsertSizeChecker.py <forward> <rev> <refGenome> <outDir> <libName> <Estimated insert size>"
            exit()
            
        logging.basicConfig(format="%(asctime)-25s%(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p")
        Configuration.instance.setOption("overwrite", "0")
        Configuration.instance.setOption("overwrite", "20")
        self.checkInsertSize(sys.argv[4], sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[5], sys.argv[6])
        Reporter.instance.createReport(sys.argv[4], small=True)
    
    def checkInsertSize(self, outputDir, forwardReads, reversedReads, refGenome, libName, insertSize):
        """
        The method checkInsertSize is the main flow of retrieving the insert size. First the sequences are mapped against
        the reference genome with bowtie. Then a histogram is created of the insert sizes. After this a R script is executed
        to draw the histogram.
        """
        self.libName = libName
        outDir = outputDir + "insertSizeEstimation/"
        if not os.path.isdir(outDir):
            os.makedirs(outDir)
        samFile = Mappers.Bowtie(outputDir, forward=forwardReads, reversed=reversedReads, refGenome=refGenome,insertSize=insertSize).execute()
        insertSizeHisto = outDir + "insertSize.histo"
        [self.mean,self.stdev, self.maxInsertSize] = InsertSizeCalculator().__getInsertSize__(samFile, insertSizeHisto)
        logging.info("mean insert size: " + str(self.mean))
        logging.info("stdev insert size: " + str(self.stdev))
        self.insertSizeHistoPlot = InsertSizePlotter(outDir, insertSizeHisto = insertSizeHisto).execute()
    
    def getLaTeXReport(self):
        txt = "\\section{Insert size estimation of "+self.libName.replace("_"," ")+"}\n"
        
        table = LaTeX.ltxTable(2)
        table.addRow(["Estimated insert size: ", str(int(round(self.mean)))])
        table.addRow(["Estimated standard deviation: ", str(int(round(self.stdev)))])
        txt = txt + table.getText()
        
        image = LaTeX.ltxImage(self.insertSizeHistoPlot)
        txt = txt + image.getText()
        
        return txt
        
    def getHtmlReport(self, reportDir):
        html = "<div class=\"insert_size\">"
        html = html + "<h1 id=\"insrtsize\">Insert size estimation of "+self.lib.libName+"</h1>"
        html = html + "<table><tr><td>Insert size: </td><td>"+str(int(round(self.mean)))+"</td></tr>"
        html = html + "<tr><td>Stdev: </td><td>"+str(int(round(self.stdev)))+"</td></tr></table>"
        html = html + "<img src=\""+os.path.relpath(self.lib.insertSizePlot, reportDir)+"\" />"
        html = html + "</div>"
        return html
    
class InsertSizePlotter(Command.Command):
    """
    This class executes the R script plot_insertsizes for plotting the insert sizes of a given histogram.
    """
    def setCommand(self):
        self.outputFile = self.outputDir + "insertSize.png"

        self.addArg("Rscript") 
        self.addArg(os.path.dirname(__file__) + "/plot_insertsizes.R")
        self.addArg(self.insertSizeHisto)
        self.addArg(self.outputFile)       

class InsertSizeCalculator(object):    
    """
    The method InsertSizeCalculator calculates the insert size for all perfect mapped reads. 
    Script edited by Jetse but mainly from: 
    http://allaboutbioinfo.blogspot.nl/2012/04/estimating-paired-end-read-insert.html
    """
    objmrl=re.compile('([0-9]+)M$');
    
    def getmeanval(self, dic,maxbound=-1):
        nsum=0;  n=0;
        for (k,v) in dic.items():
            if maxbound!=-1 and k>maxbound:
                continue;
            nsum=nsum+k*v;
            n=n+v;
        meanv=nsum*1.0/n;
        nsum=0; n=0;
        for (k,v) in dic.items():
            if maxbound!=-1 and k>maxbound:
                continue;
            nsum=nsum+(k-meanv)*(k-meanv)*v;
            n=n+v;
        varv=math.sqrt(nsum*1.0/(n-1));
        return (meanv,varv);
     
    
    
    def __getInsertSize__(self, samFile, histoFile): 
        self.plrdspan={};
        nline=0;
        for lines in open(samFile):
            field=lines.strip().split();
            nline=nline+1;
            if nline%1000000==0:
                logging.debug(str(nline/1000000)+'M...');
            if len(field)<12:
                continue;
            try:
                mrl=InsertSizeCalculator.objmrl.match(field[5]);
                if mrl==None: # ignore non-perfect reads
                    continue;
                if field[6]!='=':
                    continue;
                dist=int(field[8]);
                if dist<=0: # ignore neg dist
                    continue;
                if dist in self.plrdspan.keys():
                    self.plrdspan[dist]=self.plrdspan[dist]+1;
                else:
                    self.plrdspan[dist]=1;
            except ValueError:
                continue;
         
        if len(self.plrdspan)==0:
            logging.error('No qualified paired-end reads found. Are they single-end reads?');
        else:
            with open(histoFile, "w") as histoWriter:
                for k in sorted(self.plrdspan.keys()):
                    histoWriter.write(str(k) + "\t" + str(self.plrdspan[k]) + "\n")
            
            maxv=max(self.plrdspan,key=self.plrdspan.get);
            spanval=self.getmeanval(self.plrdspan,maxbound=maxv*3)
            return[spanval[0],spanval[1],maxv]
        
if __name__ == '__main__':
    InsertSizeChecker().main()