from Bio import SeqIO
from qualityControl import Exceptions
from collections import Counter

def fastqControl(forwardFastq, reversedFastq=None, dna=True):
    print("Checking: " + forwardFastq)
    forwardRecords = _fastqContentControl(forwardFastq, dna)
    print("Correct! File contains {} reads".format(forwardRecords))
    if reversedFastq != None:
        print("Checking: " + reversedFastq)
        reversedRecords = _fastqContentControl(reversedFastq, dna)
        print("Correct! File contains {} reads".format(reversedRecords))
        if forwardRecords != reversedRecords:
            raise Exceptions.FileFormatException("Forward file has not the same number of sequences as the reversed file when comparing " + forwardFastq + " with " + reversedFastq)
        print("Paired end is correct!")
def _fastqContentControl(fastqFile, dna=True):
    noOfRecords = 0
    try:
        for record in SeqIO.parse(open(fastqFile), "fastq"):
            lcs = Counter(record.seq.upper())
            if dna==True and sum([lcs["A"],lcs["T"],lcs["C"],lcs["G"],lcs["N"]])!=len(record.seq):
                raise Exceptions.FileFormatException("Illegal character found in " + fastqFile + " with sequence id: " + record.id)
            if max(record.letter_annotations["phred_quality"]) > 72 or min(record.letter_annotations["phred_quality"]) < 0:
                raise Exceptions.FileFormatException("Invalid quality score in " + fastqFile + " with sequence id: " + record.id)
            noOfRecords += 1
    except ValueError as e:
        raise Exceptions.FileFormatException(e)
    if noOfRecords == 0:
        raise Exceptions.FileFormatException(fastqFile + " contains no sequences...")
    return noOfRecords
