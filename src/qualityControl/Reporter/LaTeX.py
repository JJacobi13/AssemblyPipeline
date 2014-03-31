class ltxTable():
    
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
            self.text = self.text + " & " + columns[i]
        self.text = self.text + "\\\\\n"
        self.text = self.text + "\\hline\n"
        
    def getText(self):
        self.text = self.text + "\\end{tabular}\n"
        return self.text
    
class ltxImage():
    
    def __init__(self, img):
        self.img = img
        
    def getText(self):
        txt = "\\begin{figure}[h!]\n"
        txt = txt + "\\includegraphics[width=10cm]{" + self.img + "}\n"
        txt = txt + "\\end{figure}\n"
        return txt