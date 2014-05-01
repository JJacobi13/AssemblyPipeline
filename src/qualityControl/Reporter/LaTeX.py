from qualityControl import Reporter

class ltxTable():
    """
    Use this class when creating a table in a latex report. All borders and the maximum width is regulated within this class. It returns the table in LaTeX format.
    """
    def __init__(self, columns):
        self.columns = columns
        self.text = "\\begin{tabular}{| "
        for _ in range(columns):
            self.text = self.text + "l |"
        self.text = self.text + "}\n"
        self.text = self.text + "\\hline\n"
    
    def addRow(self, columns):
        if len(columns) != self.columns or len(columns) == 0:
            raise Exception("Not the same amound of columns specified")
        
        self.text = self.text + columns[0]
        for i in range(1,len(columns)):
            self.text = self.text + " & " + columns[i].replace("_"," ")
        self.text = self.text + "\\\\\n"
        self.text = self.text + "\\hline\n"
        
    def getText(self):
        self.text = self.text + "\\end{tabular}\\\\\n"
        return self.text
    
class ltxImage():
    """
    Use this class when drawing an image in a LaTeX file, it regulates the with and position of the image.
    """
    def __init__(self, img):
        self.img = img
        
    def getText(self):
        txt = "\\begin{figure}[h!]\n"
        txt = txt + "\\includegraphics[width=10cm]{" + self.img + "}\n"
        txt = txt + "\\end{figure}\n"
        return txt
    
class ltxPart():
    """
    Create this class when you need to create a new part (for instance the supplementary materials)
    """
    def __init__(self, partName):
        self.partName = partName
        Reporter.Reporter.instance.objects.append(self)
        
    def getLaTeXReport(self):
        return("\\part*{"+self.partName+"}\n\\addcontentsline{toc}{part}{"+self.partName+"}")

class ltxSection():
    """
    This class creates a new section on a new page.
    """
    def __init__(self, sectionName):
        self.txt = "\\clearpage\n"
        self.txt += "\\section{"+sectionName.replace("_"," ")+"}\n"
        Reporter.Reporter.instance.objects.append(self)
        
    def getLaTeXReport(self):
        return self.txt
        
        
        
        