from commandLineCommands import Command
from utils import FastqUtils
import os

class AllpathsAssembler(object):
    """
    The AllpathsAssembler executes an assembly with Allpaths
    """
    def doAssembly(self, outputDir, pool):
        """
        This method is the main flow of the allpaths assembly: Create the output directories, create a groups file, create a libs file, prepare the input and execute allpaths.
        """
        if not os.path.isdir(outputDir+ "assembly/"):
            os.makedirs(outputDir+ "assembly/")
        if not os.path.isdir(outputDir+ "assembly/assembly/"):
            os.makedirs(outputDir+ "assembly/assembly")
        groupsCsv = self.createInGroupsCsv(pool, outputDir + "assembly/assembly/")
        libsCsv = self.createInLibsCsv(pool, outputDir + "assembly/assembly/")
        AllpathsInputPreparator(outputDir + "assembly/assembly/", groupsCsv=groupsCsv, libsCsv=libsCsv, forward=pool.libs[0].forward).execute()
        return AllpathsExecutor(outputDir).execute()
     
    def createInLibsCsv(self, pool, outDir):
        """
        The method createInLibsCsv creates a csv file with all libraries as allpaths wants its input.
        """
        csvFile = outDir + "/libs.csv"
        with open(csvFile, "w") as csvWriter:
            csvWriter.write("library_name, project_name, organism_name, type, paired, frag_size, frag_stddev, insert_size, insert_stddev, read_orientation, genomic_start, genomic_end\n")
            for lib in pool.libs:
                csvWriter.write(lib.libName + ", ")
                csvWriter.write("assembly, ")
                csvWriter.write("assembly, ")
                if lib.type == "mp":
                    csvWriter.write("jumping, ")
                    csvWriter.write("1, ")
                    csvWriter.write(", , ")
                    csvWriter.write(str(int(lib.insertSize)) + ", ")
                    csvWriter.write(str(int(int(lib.insertSize)*0.2)) + ", ")
                    csvWriter.write("outward, ")
                elif lib.type == "pe":
                    csvWriter.write("fragment, ")
                    csvWriter.write("1, ")
                    csvWriter.write(str(int(lib.readlen)) + ", ")
                    csvWriter.write(str(int(int(lib.readlen)*0.2)) + ", ")
                    csvWriter.write(", , ")
                    csvWriter.write("inward, ")
                elif lib.type == "u":
                    csvWriter.write("long, ")
                    csvWriter.write("0, ")
                    csvWriter.write(", , , , , ")
                csvWriter.write("0, 0\n")
        return csvFile
    
    def createInGroupsCsv(self, pool, outputDir):
        """
        The method createInGroupsCsv creates a csv file with all paths to the libraries in it. When the reads are mated, the last "1" is replaced by a ? for the allpaths regex input.
        """
        csvFile = outputDir + "/groups.csv"
        with open(csvFile, "w") as csvWriter:
            i = 0
            csvWriter.write("group_name, library_name, file_name\n")
            for lib in pool.libs:
                if lib.reversed != None:
                    csvWriter.write(str(i) + ", " + lib.libName + ", " + rreplace(lib.forward, "1", "?",1) + "\n")
                    i += 1
                else:
                    csvWriter.write(str(i) + ", " + lib.libName + ", " + lib.forward + "\n")
                    i += 1
                    
        return csvFile
def rreplace(s, old, new, occurrence):
    """
    The method rreplace replaces the last occurence(s) of a given string with a new value.
    :param s: The string to replace the last value in
    :param old: The substring to replace
    :param new: The new value for the substring
    :param occurence: The numer of times to replace the substring, beginning from the end.
    """
    li = s.rsplit(old, occurrence)
    return new.join(li)
               
class AllpathsInputPreparator(Command.Command):
    """
    The AllpathsInputPreparator executes a script given with allpaths to prepare the input of allpaths.
    """

    def setCommand(self):
        self.outputFile = self.outputDir + "jump_reads_orig.fastb"
        self.addArg("PrepareAllPathsInputs.pl")
        self.addArg("DATA_DIR=" + self.outputDir)
        self.addArg("IN_GROUPS_CSV=" + self.groupsCsv)
        self.addArg("IN_LIBS_CSV=" + self.libsCsv)
        if FastqUtils.determineQuality(self.forward) == 64:
            self.addArg("PHRED_64=True")
        self.addArg("PLOIDY=2")
        self.addArg("OVERWRITE=1")
               
class AllpathsExecutor(Command.Command):
    """
    This command executes allpaths.
    """
    def setCommand(self):
        self.outputFile = self.outputDir + "/assembly/assembly/assembly/ASSEMBLIES/assembly/final.contigs.fasta"
        self.scaffolds = self.outputDir + "/assembly/assembly/assembly/ASSEMBLIES/assembly/final.assembly.fasta"
        self.addArg("RunAllPathsLG")
        self.addArg("PRE=" + self.outputDir)
        self.addArg("REFERENCE_NAME=assembly")
        self.addArg("DATA_SUBDIR=assembly")
        self.addArg("RUN=assembly")
        self.addArg("SUBDIR=assembly")
        self.addArg("OVERWRITE=TRUE")
        self.addArg("FF_MAX_STRETCH=5") #TODO FIND OUT WHAT THIS PARAMETER MEANS!
        self.addArg("THREADS=" + Configuration.instance.getGlobalOption("maxThreads"))
        
        
if __name__ == '__main__':
    from model import Pool
    from configuration import Configuration
    import sys, logging
    
    logging.basicConfig(filename="pipeline.log", format="%(asctime)-25s%(levelname)-12s%(message)s", level=logging.DEBUG, datefmt="%m/%d/%Y %I:%M:%S %p")
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter("%(asctime)-25s%(message)s", "%I:%M:%S %p"))
    logging.getLogger().addHandler(console)
    Configuration.instance.parseIni(sys.argv[1])
    Configuration.instance.setOption("overwrite", "1")
    pool = Pool.Pool(Configuration.instance.getGlobalOption("outputDirectory"))
    pool.createLibs(Configuration.instance.libNames)
    AllpathsAssembler().doAssembly(Configuration.instance.getGlobalOption("outputDirectory"), pool)
        