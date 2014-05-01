'''
Created on Jun 21, 2013

@author: Jetse
'''
from model import Library, Collection
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

    def createLibs(self, libNames):
        for lib in libNames:
            libInstance = Library.Library(self, lib, Configuration.instance.getLibInfo(lib))
            self.libs.append(libInstance)
