import re, sys

class GenBlastParser():
    def __init__(self, filename):
        self.currentGene = None
        self.genes = []
        self.queries = {}
        with open(filename) as genblastReader:
            for line in genblastReader:
                if line.startswith("//") or line.startswith("genBlast phase"):
                    if line.startswith("//for query"):
                        line = line.replace("//","").strip()
                        self.queries[line.split("|")[2]] = [0,0]
                    continue
                if line.strip():
                    self.parseGenes(line)
        partial = 0
        full = 0
        for entry in self.queries.values():
            if entry[0] > 0:
                full += 1
            if entry[1] > 0 or entry[0] > 0:
                partial += 1
                        
        self.partial = partial / float(len(self.queries)) * 100
        self.full = full/float(len(self.queries))*100
    
    def parseGenes(self, line):
        splittedLine = line.split("|")
        if len(splittedLine) == 8:
            self.currentGene = Gene(splittedLine)
            if float(self.currentGene.score) > 0:
                if self.currentGene.cover > 70:
                    self.queries[splittedLine[2]][0] += 1
                else:
                    self.queries[splittedLine[2]][1] += 1
                self.genes.append(self.currentGene)
            else:
                self.currentGene = None
        else:
            if self.currentGene !=None:
                try:
                    Exon(line, self.currentGene)
                except StandardError:
                    #This line probaby contains None, so this gene has no hits
                    pass
    
    def createGffFile(self, fileName):
        with open(fileName, "w") as gffWriter:
            for gene in self.genes:
                gffWriter.write(gene.toGff())
        
class Gene():
    def __init__(self, geneInfo):
        self.chrom = geneInfo[3].split(":")[0]
        self.start = geneInfo[3].split(":")[1].split("..")[0]
        self.end = geneInfo[3].split(":")[1].split("..")[1]
        self.name = geneInfo[2] + "_" + geneInfo[7].split(":")[1].strip()
        self.score = geneInfo[6].split(":")[1]
        self.rank = int(geneInfo[7].split(":")[1].strip())
        self.cover = float(re.search("\((\d+\.?\d*)%\)",geneInfo[5]).group(1))
        self.exons = []
    
    def toGff(self):
        txt = "\t".join([self.chrom,"genBlastA","gene",self.start,self.end,self.score, ".", ".","name="+self.name+";"]) + "\n"
        for exon in self.exons:
            txt += exon.toGff()
        return txt
    
    def __repr__(self):
        return self.__str__()
        
    def __str__(self):
        return "Gene[name="+self.name+ ", chrom=" + self.chrom + ", start=" +self.start + ", end=" + self.end + ", score="+self.score+ ", exons=" + str(len(self.exons))+"]"
        
    
        
class Exon():
    
    def __init__(self, exonInfo, gene):
        matchObj = re.search("HSP_ID\[(\d+)\]:\(([0-9]+).([0-9]+)\);query:\(([0-9]+).([0-9]+)\); pid: (\d+\.?\d*)", exonInfo)
        if matchObj == None:
            raise StandardError("No match found...")
        self.exonIndex = matchObj.group(1)
        self.start = matchObj.group(2)
        self.end = matchObj.group(3)
        self.identity = matchObj.group(6)
        gene.exons.append(self)
        self.gene = gene
    
    def toGff(self):
        return "\t".join([self.gene.chrom,"genBlastA","exon",self.start,self.end,self.identity, ".", ".","parent="+self.gene.name+";"]) + "\n"
        
    def __repr__(self):
        return self.__str__()
        
    def __str__(self):
        return "Exon[id="+self.exonIndex+ ", start=" +self.start + ", end=" + self.end + ", identity="+self.identity+"]"
    
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "usage: <genblast file> <output gff file>"
        exit()
    gbp = GenBlastParser(sys.argv[1])
    gbp.createGffFile(sys.argv[2])