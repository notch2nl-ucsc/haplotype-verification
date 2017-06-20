#!/usr/bin/env python2.7

# This script collects original FASTQ sequence/quality pairs based on their names. Used in conjunction with the paralog
#   binning script to see whether full sequence quality is a determining factor in poorly matching reads.
#
# Average read quality (assuming illumina 1.8 phred+33) is printed to stdout
#
# Input Files:
#   -FASTQ containing all reads
#   -FASTA file containing a subset of reads from the FASTQ


referenceFastaFilename = "RESULTS_margin_align_longer_consensus/exceptionalReads.fa"
fastqFilename = "interesting_seqs.fq.txt"
outputFastqFilename = "exceptionalReads.fq"

readNames = set()
with open(referenceFastaFilename, 'r') as referenceFile:
    for line in referenceFile:
        if line.startswith('>'):
            readNames.add(line[1:].split()[0])

readQualities = dict()
readSequences = dict()
qualities = list()
with open(fastqFilename, 'r') as fastqFile:
    readName = None
    isQual = False
    isSeq = False
    for line in fastqFile:
        if isSeq==True:
            readSequences[readName] = line.strip()
            isSeq=False

        if line.startswith('@'):
            readName = line.strip()[1:].split()[0]
            # if readName in readNames:
                # print(readName)
            isSeq=True

        if isQual==True:
            if readName in readNames:
                qString = line.strip()
                # print(qString)

                for char in qString:
                    q = ord(char) - 33
                    # print(q)
                    qualities.append(q)
                readQualities[readName] = qString
            isQual=False

        if line.startswith('+'):
            isQual = True
            isSeq=False

print("Average quality: %d"%(sum(qualities)/len(qualities)))


with open(outputFastqFilename,'w') as outFile:
    for read in readQualities:
        name = read
        quality = readQualities[name]
        sequence = readSequences[name]
        outFile.write("@%s\n%s\n+\n%s\n"%(read,sequence,quality))

