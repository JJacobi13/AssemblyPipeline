# from model import Command
# import os
# 
# class FastqConverter(object):    
#     
#     def convertToFastq(self, lib):
#         if lib.format == "sff":
#             sffToFastqConverter(lib).execute()
#             
#             
# class sffToFastqConverter(Command.Command):
#     
#     def setCommand(self):
#         self.addArg("sff2fastq")
#         self.addArg("-o " + self.collection.outputDir + self.collection.rawDir)
#         
# 
#     def updateStatus(self):
#         pass
#             
#     def concatFastq(self):
#         if self.path.find(",") != -1:            
#             outFile = self.outputDir + "rawInput.fastq"
#             Program.Program.execute("cat " + self.path.replace(",", " ") + ">" + outFile, self)
#             self.path = outFile
#          
#     def convertSffTofastq(self):
#         files = self.path.split(",")
#         self.path = ""
#         for origFile in files:
#             fastq=os.path.splitext(os.path.basename(origFile))[0]+".fastq"
#             self.path=self.path + self.outputDir + fastq + ","
#             Program.Program.execute(" -o " + self.outputDir + fastq + " " + origFile, self)
#         self.path = self.path[:-1]
#         self.format = "fastq"
            