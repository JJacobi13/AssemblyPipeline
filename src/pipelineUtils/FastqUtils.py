# import sys
# from itertools import izip
def determineQuality(fastqFile):
    """
    The method determineQuality determines whether a given fastq file is phred 33 or phred 64 encoded.
    """
    lineNo = 0;
    lowest = 100;
    highest = 0;
    with open(fastqFile) as fastqReader:
        for line in fastqReader:
            lineNo += 1
            if lineNo % 4 == 0:
                for char in line:
                    lowest = min([lowest,ord(char)])
                    highest = max([highest, ord(char)])
                    
            if lineNo > 1000:
                if highest > 77:
                    return 64
                return 33

# def mergeFastqFiles(forwardReads, reversedReads, outFile):
#     revHash = []
#     i = 0
#     with open(forwardReads) as ffh, open(reversedReads) as rfh:
#         with open(outFile, "w") as writer:
#             for fr, rr in izip(ffh, rfh):
#                 writer.write(fr)
#                 revHash.append(rr)
#                 if len(revHash) == 4:
#                     i += 1
#                     if i % 1000000 == 0:
#                         print "Merged " + str(i)/1000000 + "m reads"
#                     for line in revHash:
#                         writer.write(line)
#                         revHash = []
#     
#             
# if __name__ == '__main__':
# #     mergeFastqFiles("D:/Parasponia/miseq_PE_400_subset_1.fastq", "D:/Parasponia/miseq_PE_400_subset_2.fastq", "D:/Parasponia/merged.fastq")
#     if sys.argv[1] == "qual":
#         print determineQuality(sys.argv[2])
#     elif sys.argv[1] == "merge":
#         mergeFastqFiles(sys.argv[2],sys.argv[3], sys.argv[4])