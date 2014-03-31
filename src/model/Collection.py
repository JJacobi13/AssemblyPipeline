import os, stat

class Collection(object):
    """
    The collection is an abstract object with methods used in both library objects and pool objects.
    """
    
    def setup(self, outputDir):
        """
        The setup method creates the output dir of the collection. This method also creates a job file with all commands of this collection.
        
        :param outputDir: The directory this collection writes its output to
        :fileType pool: str (filepath)
        """
        self.outputDir = outputDir
        if not os.path.isdir(self.outputDir):
            os.makedirs(self.outputDir)
             
        self.jobFile = self.outputDir + "/" + "job.sh"
        
        with open(self.jobFile, "w") as jobWriter:
                jobWriter.write("#!/bin/bash" + "\n")
                
        os.chmod(self.jobFile, stat.S_IRWXU)
            
    def addCommand(self, command):
        """
        This method is always called when a program is executed. It adds the command to the job file
        
        :param command: the command to add to the shell file
        :fileType command: str
        """
        with open(self.jobFile, "a") as jobWriter:
            jobWriter.write(command + "\n")