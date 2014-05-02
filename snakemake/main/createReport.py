import glob
from qualityControl.Reporter import Reporter
from qualityControl import FastqSmallReport
libs = {}
for fileName in glob.glob("*.fastq"):
    lib = fileName.split("_")[0].split(".")[0]
    if lib not in libs:
        libs[lib] = []
    libs[lib].append(fileName)

for key, value in libs.items():
    value.sort(key=lambda item: (len(item), item))
    i = 0
    while i < len(value):
        forwardFastq = value[i]
        reversedFastq = None  
        try:
            if  len(forwardFastq) == len(value[i+1]):
                reversedFastq = value[i+1]
                i += 1
        except IndexError:
            pass
        i += 1
        if reversedFastq == None:
            print("Creating report for " + forwardFastq)
        else:
            print("Creating report for " + forwardFastq + " and " + reversedFastq)
        FastqSmallReport.FastqSmallReport().createSmallReport(forwardFastq, reversedFastq)
        Reporter.Reporter().createReport("report/", small)
        
