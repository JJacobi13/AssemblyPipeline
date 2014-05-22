'''
Created on May 16, 2014

@author: Jetse
'''

#First with .fastq.gz, because .gz gives maximum recursion depth exceeded... 
#TODO: find better way to extract all gz files instead of .fastq.gz only...
rule unzip:
    input: "{file}.fastq.gz"
    output: "{file}.fastq"
    shell: "gunzip {input[0]}"