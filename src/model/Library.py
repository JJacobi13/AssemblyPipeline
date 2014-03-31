from model import Collection

class Library(Collection.Collection):
    """
    The library class represents a Library.
    
    """
    
    def __init__(self, pool, libName, libInfo=None):
        """
        The constructor of the Sample object sets the given parameters as instance variables and creates the output directory if not exists
        
        :param pool: the pool this sample is a part of
        :fileType pool: :py:class:`Pool.Pool`
        :param libName: The name of this sample
        :fileType libName: Str
        """
        
        self.pool = pool
        self.libName = libName
        if libInfo != None:
            self.rawForward = libInfo.getOption("forward")
            self.rawReversed = libInfo.getOption("reversed")
            self.forward = libInfo.getOption("forward")
            self.reversed = libInfo.getOption("reversed")
            self.format = libInfo.getOption("format")
            self.type = libInfo.getOption("type")
            self.insertSize=libInfo.getOption("insertSize")
            self.stdev=libInfo.getOption("stdev")
            self.readlen=libInfo.getOption("readlen")
            self.sequencingPlatform=libInfo.getOption("sequencingPlatform")
        self.status = "raw"  
        self.setup(pool.outputDir + libName + "/")   

    def __str__(self):
        return "Sample[Lib name: " + self.libName + ", forward: " + str(self.forward) + ", reversed: " + str(self.reversed) + " ]"

    def __repr__(self):
        return self.__str__()