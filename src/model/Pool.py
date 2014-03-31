'''
Created on Jun 21, 2013

@author: Jetse
'''
import logging

from model import Library, Collection
from time import localtime, strftime
from configuration import Configuration

class Pool(Collection.Collection):
    """The Pool object represents a pool of samples where programs can be executed on
    
    """

    def __init__(self, outputDir):
        """
        The constructor of the pool creates a configuration object. ONLY GIVE NO CONFIGFILE FOR TESTING OR EXECUTING CLASSES INDEPENDEND OF THE PIPELINE!!!
        """
        
        self.libs = []
        self.setup(outputDir + "/")
        
        self.configureLogging()
    
    def configureLogging(self):
        """Setting up the variables for the logging module.
        
        """
        logging.basicConfig(filename=self.outputDir + strftime("%d-%m-%y_%H-%M", localtime()) +".log", format="%(asctime)-25s%(levelname)-12s%(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p")
        #logging.basicConfig(filename="C:/Users/Jetse/workspace/assembly/assembly.log", format="%(asctime)-22s%(levelname)-12s%(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p")
        
        console=logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(logging.Formatter("%(asctime)-25s%(message)s", "%I:%M:%S %p"))
        logging.getLogger().addHandler(console)
        

    def createLibs(self, libNames):
        for lib in libNames:
            libInstance = Library.Library(self, lib, Configuration.instance.getLibInfo(lib))
            self.libs.append(libInstance)
