#!/usr/bin/env python2.7

# This script creates binned file dumps, using MarkovAligner to bin and score FASTA reads:
# Each haplotype has a set of files with "TiesExcluded" or "TiesSplit" suffixes, which correspond to whether ambiguous
# reads were thrown out. Given a haplotype+suffix I have the following files:
#       .fa    -    a list of the unaligned reads that were estimated to match that haplotype with their bit score
#                   (from -3 to 0) tagged onto their name
#       .txt   -    a feature table for all the aligned reads assigned to the haplotype (with the first row being the
#                   reference haplotype, see far right columns for names/scores)

from MarkovAligner import *
from collections import defaultdict
from matplotlib import pyplot as plot
import vcf
import pandas as pd
import pysam
import os

# samFileName = "alignedReadsFrankANDYuliaANDMeghan.sam"
samFileName = "margin_align_longer_consensus.sam"
# samFileName = "alignedReadsMeghan.sam"

outFolder = "RESULTS_"+samFileName.split('.')[0]

if not os.path.exists(outFolder):
    os.makedirs(outFolder)

readsByName = dict()

with open(samFileName) as samfile:
    for line in samfile:
        if not line.startswith('@'):
            line = line.split('\t')
            readID = line[0]
            readSequence = line[9]

            # print(readID)

            readsByName[readID] = readSequence


#----------------------------------------------------------------------------------------------------------------------#
# FROM IAN FIDDES

reads = list(pysam.Samfile(samFileName))

recs = list(vcf.Reader(open('variants_longer.vcf')))
recs = [x for x in recs if x.POS <= 1100]
intervals = []

for r in recs:
    start, stop = r._compute_coordinates_for_indel()
    intervals.append([start - 1, stop, r])  # convert to 0-based

# how many distinct feature vectors are there?
feature_vectors = defaultdict(list)  # why defaultdict? <- initializes all values as an empty list
for r in recs:
    for sample in r.samples[:-1]:  # exclude consensus
        gt = sample.gt_bases.split('/')[0]  # remove duplicate SNP (homozygous)
        feature_vectors[sample.sample].append(gt)  # store sample name:[SNP1,SNP2...SNPn] in dictionary

feature_df = pd.DataFrame.from_dict(feature_vectors).T
feature_df = feature_df.drop_duplicates()

# print(feature_df)

results = []
for a in reads:  #nanopore MarginAlign-ed reads
    r = [a.qname]
    m = {y: x for x, y in a.get_aligned_pairs()}  #a list of aligned read (query) and reference positions
    for start, stop, rec in intervals:  #only collect the SNPs/features from these alignments
        start = m.get(start, None)
        stop = m.get(stop, None)
        if start is not None and stop is not None:
            seq = a.seq[start:stop]
        else:
            seq = None
        r.append(seq)  #store the features as a list
    results.append(r)

#----------------------------------------------------------------------------------------------------------------------#


names = feature_df.index.tolist()[:]        # names of each haplotype (from Alex)
featureArray = feature_df.as_matrix()[:,:]  # set of features for the proposed haplotypes

# ambiguous reads excluded
variantCounts_noTies = {k: 0 for k in names}    # dataset-wide counts of each alignment
variantBins_noTies = {k: list() for k in names} # contains (readName,score) for each aligned read

# ambiguous reads split evenly
variantCounts_all = {k: 0 for k in names}
variantBins_all = {k: list() for k in names}


model = MarkovAligner()
model.buildModel(featureArray,names)    #generate a model based on the proposed haplotypes


allFeatureVectors = [x[1] for x in results]
nVectors = len(allFeatureVectors)
# nTraining = nVectors*2/3    #unused at the moment

print("Total Feature Vectors: %d"%nVectors)
# print("Training Vectors: %d"%nTraining)

# trainingFeatureVectors = allFeatureVectors[:nTraining]
# testingFeatureVectors = allFeatureVectors[nTraining:]


for entry in results[:]:
    readName = entry[0]
    readFeatureVector = entry[1:]

    # find the most probable haplotype(s)
    paths,maxPosteriorProb = model.maximumProbabilityPath(readFeatureVector)

    for path in paths:
        haplotypeMatchName = path

        # track alignment counts for each haplotype (splitting ties)
        variantCounts_all[haplotypeMatchName] += 1.0/len(paths)
        variantBins_all[haplotypeMatchName].append((readName,maxPosteriorProb))

        if len(paths) == 1: # (ties not allowed)
            variantCounts_noTies[haplotypeMatchName] += 1.0
            variantBins_noTies[haplotypeMatchName].append((readName,maxPosteriorProb))
            # alignmentScores_noTies.append((readName,maxPosteriorProb))


print("\nVariant Total Counts (split ties):")

total = 0
for key in sorted(variantCounts_all):
    count = variantCounts_all[key]
    total += count
    print("%s\t%.2f"%(key,count))
print("TOTAL\t%.2f"%total)


print("\nVariant Total Counts (ties excluded):")

total = 0
for key in sorted(variantCounts_noTies):
    count = variantCounts_noTies[key]
    total += count
    print("%s\t%d"%(key,count))
print("TOTAL\t%d"%total)


binnedScoresAll = dict()

#create dump files
for binName in variantBins_all:
    binnedReads = list(variantBins_all[binName])
    binScores = list()

    with open(outFolder+'/'+binName+"_alignedReadsTiesSplit.txt", 'w') as binFile,\
         open(outFolder+'/'+binName+"_alignedReadsTiesSplit.fa", 'w') as fastaFile:

        binFile.write('\t'.join(feature_df.loc[binName])+"\tASSEMBLED %s\n"%binName)

        if len(binnedReads) != 0:
            sortedResults = sorted(binnedReads,key=lambda x: x[1],reverse=True)

            for entry in sortedResults:
                name = entry[0]
                score = entry[1]

                binScores.append(2**score)

                featureVector = [x[1:] for x in results if x[0] == entry[0]][0]  #inefficient search
                featureVector = '\t'.join([x if x!=None else 'N' for x in featureVector])

                fullRead = readsByName[entry[0]]

                binFile.write(featureVector +'\t'+ '\t'.join(map(str,list(entry))) +'\n')

                fastaFile.write(">%s_%.4f\n%s\n"%(name,score,fullRead))

    binnedScoresAll[binName] = binScores


binnedScoresNoTies = dict()

for binName in variantBins_noTies:
    binnedReads = list(variantBins_noTies[binName])
    binScores = list()

    # print(binnedReads)

    with open(outFolder+'/'+binName+"_alignedReadsTiesExcluded.txt", 'w') as binFile,\
         open(outFolder+'/'+binName+"_alignedReadsTiesExcluded.fa", 'w') as fastaFile:

        binFile.write('\t'.join(feature_df.loc[binName])+"\tASSEMBLED %s\n"%binName)

        if len(binnedReads) != 0:
            sortedResults = sorted(binnedReads,key=lambda x: x[1],reverse=True)

            binScores.append(sortedResults[0][1])

            for entry in sortedResults:
                name = entry[0]
                score = entry[1]


                featureVector = [x[1:] for x in results if x[0] == entry[0]][0]  #inefficient search
                featureVector = '\t'.join([x if x!=None else 'N' for x in featureVector])

                fullRead = readsByName[entry[0]]

                binFile.write(featureVector +'\t'+ '\t'.join(map(str,list(entry))) +'\n')

                fastaFile.write(">%s_%.4f\n%s\n"%(name,score,fullRead))

        else:
            pass

    binnedScoresNoTies[binName] = binScores


readsSortedByMissingValues = list()

#Create a read dump sorted by a crude measurement of quality (number of missing positions in the feature vector)
for entry in results[:]:
    read = [x if x!=None else '*' for x in entry[1:]]
    quality = len([x for x in read if x!=None])

    readsSortedByMissingValues.append((read,quality))

readsSortedByMissingValues.sort(key=lambda x: x[1])


fig,axes = plot.subplots(len(names),1,sharex=True,sharey=True)

# plot the score distributions as log histograms for each haplotype's bin
for h,bin in enumerate(sorted(binnedScoresAll)):
    x = binnedScoresAll[bin]
    axes[h].hist(x,25,normed=0,histtype='bar',color=[.118, .565, 1.00],linewidth=0.1,range=[0.125,1.0],log=True)
    axes[h].set_ylabel(bin[-3:])

    if h == len(axes)-1:
        axes[h].set_xlabel("Match Likelihood")


plot.show()


