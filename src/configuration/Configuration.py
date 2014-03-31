"""
The configuration module regulates all configuration settings which are in the ini file.
This is done by 2 objects, :class:`Config` and :class:`Section`.
The :class:`Config` regulates the parsing of the ini file and the communication with other objects
The :class:`Section` represents a section of the ini file
"""
import os

class Config(object):
    """
    The Config class regulates everything of the configuration of the programs.
    This information is read from the iniFile config.ini which is located in this directory.
    """
    def __init__(self):
        self.sections = {}
        self.libNames = []
        
    def parseIni(self, iniFile):
        name="default"
        try:
            with open(iniFile) as ini:
                content = ini.readlines()
            for line in content:
                line = line.rstrip()
                if line.startswith("#"):
                    continue
                if "[" in line and "]" in line:
                    name = line.replace("[","").replace("]","")
                    self.sections[name] = Section(name)
                elif line:
                    splitLine = line.partition("=")
                    if name=="library" and splitLine[0]=="libraryname":
                        self.sections[splitLine[2]] = self.sections.pop(name)
                        name = splitLine[2]
                        self.libNames.append(name)
                    self.sections[name].setOption(splitLine[0], splitLine[2])
        except(IOError):
            raise IOError("can not read configuration file at " + os.path.abspath(iniFile))
        
    def getGlobalOption(self, value):
        return self.sections["global"].getOption(value)
    
    def setOption(self, index, value):
        if (self.sections.has_key("global") == False):
            self.sections["global"] = Section("global")
        self.sections["global"].setOption(index, value)
    
    def getLibInfo(self, libName):
        return self.sections[libName]
    
    
    def getPath(self, program):
        return self.sections["path"].GetOption(program)
                
class Section:
    """
    The class section represents a section from the ini file.
    This section contains a header name and a dictionary of key-value pairs.
    The value is an empty string if this option has no value 
    """
    def __init__(self, name):
        self.name = name
        self.options = {}
        
    def getOption(self, name):
        """
        The method getOption is a getter for retrieving a specific option      
        """
        if name in self.options and self.options[name].strip():
            return self.options[name]
        return None
    
    def setOption(self, name, value):
        """
        The method setOption is a setter for setting a new option, overwrites when exist.
        """
        self.options[name] = value
        
    def getAllOptions(self):
        """
        Getter for retrieving the full dictionary in once.
        """
        return self.options
    
instance = Config()