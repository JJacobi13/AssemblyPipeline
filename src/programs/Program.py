# from subprocess import Popen, PIPE
# import logging

class Program(object):
    """The Program is a class which represents an executable program.

    """
    @staticmethod
    def execute(command, collection):
        collection.addCommand(command)
#         logging.info("starting " + programName)
#         logging.debug("command: " + command)
# 
#        
#         error,output = Popen(command, shell=True, stdout=PIPE, stderr=PIPE).communicate()
#         if(len(output) > 0):
#             logging.info(output)
#         if(len(error) > 0):
#             logging.info(error)
#         logging.info("finished " + programName)
# 
#     @staticmethod
#     def getProgramArguments(programName):
#         """The method getProgramArguments asks the configuration object for all arguments and parses these to one space separated string.
#         When the program needs the arguments in a different way, this method can be overrided.
#         
#         :param programName: The name of the program where to retrieve the options from
#         :type programName: str.
#         :returns: The program arguments as a <space> separated string.
#         
#         """
#         args = config.getProgramOptions(programName)
#         argsCmd = " "
#         for arg in args:
#             argsCmd = argsCmd + arg + " "
#             if args[arg]:
#                 argsCmd = argsCmd + args[arg] + " "
#         return argsCmd
