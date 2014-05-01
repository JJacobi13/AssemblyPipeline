import glob

rule createSmallReport:
    input: "{samples}.fastq"
    output: "{samples}.fastq.summary"
    run:
        FastqSmallReport.FastqSmallReport().createSmallReport(input[0])
        
rule report:
    input: glob.glob("*.fasta")
    output: "report.txt"
    run:
        fastqSummaries = []
        for file in glob.glob("*.fastq"):
            FastqSmallReport.FastqSmallReport().createSmallReport(file)
            
            
        print(fastqSummaries)