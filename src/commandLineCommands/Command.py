from abc import ABCMeta, abstractmethod
import logging, subprocess, os
from configuration import Configuration

class Command(object):
    """
    The command object is the superclass of all executable commands. 
    This objects contains methods for easy creation and execution of command-line tools.
    """
    __metaclass__ = ABCMeta
    
    def __init__(self, outputDir, **kwargs):   
        """
        The constructor of the command class sets the collection as an instance variable. The constructor also creates an empty command
        which is filled with the setCommand method of the child class. After setting the command the status of the collection will be
        updated and the command will be added to the job file of the collection.
        """    
        self.outputDir = outputDir   
        if os.path.exists(outputDir) == False:
            os.mkdir(outputDir)
            
        #Set all other arguments as instance variables
        self.__dict__.update(kwargs)
        del kwargs
        self.__dict__.update(locals())
        del self.self
        
        self.cmd = ""
        self.checkInput()
        self.setCommand()
        self.updateStatus()
    
    def addArg(self, cmdToAdd):
        """
        addArg is a method for fast adding arguments to a program, without having trouble with cloting arguments.
        """
        self.cmd = self.cmd + cmdToAdd + " "
    

    def checkInput(self):
        """
        The method checkInput has to check whether the command can be executed, else a previous program has to be executed
        """
        pass
    
    @abstractmethod
    def setCommand(self):
        """
        The method setCommand retrieves a command to execute on the commandline.
        """
        pass
    
    def updateStatus(self):
        """
        The method updateStatus updates the collection to its new status.
        """
        pass
    
    def execute(self):
        from qualityControl.Reporter import Reporter
        Reporter.instance.objects.append(self)
        if int(Configuration.instance.getGlobalOption("overwrite")) == 0:
            if os.path.isfile(self.outputFile):
                logging.info("skipping command " + self.cmd.split(" ")[0] + ", outputfile already exists.")
                logging.debug("command: " + self.cmd)
                return self.outputFile
        logging.info("starting command " + self.cmd.split(" ")[0])
        logging.debug("command: " + self.cmd)

        p = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in iter(p.stdout.readline, ''):
            if(len(line) > 0):
                logging.debug(line.rstrip('\n'))
        p.wait()
        logging.info("finished " + self.cmd.split(" ")[0])
        
        return self.outputFile