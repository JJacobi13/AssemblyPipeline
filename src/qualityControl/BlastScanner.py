#results = "blast_out_1.tab"
import subprocess,  os, logging, copy
import asciitree
from configuration import Configuration

class BlastScanner(object):
    
    def getLaTeXReport(self):
        txt = "\\section{Tree of found homolog sequences}\n"
        txt = txt + "\\begin{verbatim}\n"
        txt = txt + asciitree.draw_tree(self.tree["1"])
        txt = txt + "\\end{verbatim}\n"
        return txt
    
    def scan_results(self, inputFile):
        """
        
        """
        logging.info("reading blast output")
        organismsGi = {} 
        with open(inputFile) as blastReader:
            for line in blastReader:
                info = line.split()
                if info[4] in organismsGi:
                    organismsGi[info[4]] = organismsGi[info[4]] + 1
                else:
                    organismsGi[info[4]] = 1
        
        logging.info("convert gi to tax")
        gi2taxId = Gi2taxIdConverter(Configuration.instance.getGlobalOption("giToTax"))
        gi2taxIdDict = gi2taxId.convert(organismsGi.keys())
        organismsTax = {}
        for gi in organismsGi.keys():
            if gi in gi2taxIdDict:
                organismsTax[gi2taxIdDict[gi]] = organismsGi[gi] 
                del organismsGi[gi]
        del organismsGi
        
        logging.info("creating tree")
        self.tree = {}
        self.tree["1"] = TaxonomyNode("1", 0)
        
        self.fullTax = self.parseTax()
        for taxId in organismsTax.keys():
            self.appendToTree(taxId, organismsTax[taxId])
        
        self.tree["1"].updateName(self.parseNames())
        with open(os.path.splitext(inputFile)[0] + "TreeUnpruned.txt", "w") as writer:
            writer.write(asciitree.draw_tree(self.tree["1"]))
            
        logging.info("pruning tree")
        
        self.tree["1"].pruneCounts()
        self.tree["1"].pruneParents()
        logging.info("drawing tree")
        
        
        with open(os.path.splitext(inputFile)[0] + "Tree.txt", "w") as writer:
            writer.write(asciitree.draw_tree(self.tree["1"]))
         
    def appendToTree(self, taxId, counts):
        node = TaxonomyNode(taxId, counts)
        self.tree[taxId] = node
        while True:
            if taxId in self.fullTax:
                parent = self.fullTax[taxId]
                node = self.tree[taxId]
                if parent in self.tree:
                    self.tree[parent].children.append(node)
                    node.parent = self.tree[parent]
                    self.tree[parent].updateCounts(counts)
                    return
                else:
                    self.tree[parent] = TaxonomyNode(taxId, counts)
                    node.parent = self.tree[parent]
                    self.tree[parent].children.append(node)
                taxId = parent
            else:
                self.tree["1"].children.append(TaxonomyNode(taxId))
                
                return
    
    def parseTax(self):
        taxonomy = {}
        with open(Configuration.instance.getGlobalOption("taxNodes")) as treeReader:
            for line in treeReader:
                info = line.split("|")
                taxonomy[info[0].strip()] = info[1].strip()
        return taxonomy
    
    def parseNames(self):
        names = {}
        with open(Configuration.instance.getGlobalOption("taxNames")) as treeReader:
            for line in treeReader:
                if "scientific name" in line:
                    info = line.split("|")
                    names[info[0].strip()] = info[1].strip()
        return names
        

class TaxonomyNode(object):
    
    def __init__(self, taxId, counts):
        self.parent = None
        self.taxId = taxId
        self.children = []
        self.counts = counts
        
    def updateName(self, names):
        self.name = names[self.taxId]
        for child in self.children:
            child.updateName(names)
        
    def updateCounts(self, counts):
        self.counts = self.counts + counts
        if self.parent != None:
            self.parent.updateCounts(counts)
        
    def __repr__(self):
        return self.__str__()
    
    def __str__( self ):
        if hasattr(self, "name"):
            return self.name + "(" + str(self.counts) + ")"
        else:
            return self.taxId + "(" + str(self.counts) + ")"
        
    def toString(self, depth=""):
        string = "\\-"
        string = string + self.taxId + "(" + str(self.counts) + ")\n"
        for child in self.children:
            string = string + child.toString("")
        return string
        
    def pruneParents(self):
        if self.parent !=None:
            if self.parent.parent !=None:
                if self.parent.parent.parent !=None:
                    if self.counts == self.parent.counts:
                        for index, item in enumerate(self.parent.parent.children):
                            if item == self.parent:
                                self.parent.parent.children[index] = self
                        self.parent = self.parent.parent

        for child in self.children:
            child.pruneParents()
    
    def pruneCounts(self):
        for child in copy.copy(self.children):
            child.pruneCounts() 
            
        if self.counts < 5:
            self.parent.children.remove(self)

                
    def totalNoOfElements(self):
        counter = 1
        for child in self.children:
            counter = counter + child.totalNoOfElements()
        return counter

class Gi2taxIdConverter(object):
    
    def __init__(self, inputFile):
        self.inputFile = inputFile
        
    def convert(self, ids):
        gi2tax = {}
        foundVals = []
        no = 0
        cmd = "zgrep \"^\(" + "\|".join(ids) + "\)[[:space:]]\" " + self.inputFile
        
#         cmd = "ag --no-numbers \"^(" + "|".join(ids) + ")\\s" + "\" " + self.inputFile

        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in iter(p.stdout.readline, ''):
            info = line.split()
            gi2tax.update({info[0]:info[1]})
            no = no + 1
            foundVals.append(info[0])
        p.wait()    
        self.UpdatedGi = [x for x in ids if x not in set(foundVals)]
        return gi2tax
    
# BlastScanner().scan_results("/home/jaco001/fungiAssembly/data/assembly9/contaminationCheck/blast_out.csv")