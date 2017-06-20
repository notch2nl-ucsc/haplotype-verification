#!/usr/bin/env python2.7

# This script creates binned file dumps, using MarkovAligner to bin and score FASTA reads:
# Each haplotype has a set of files with "TiesExcluded" or "TiesSplit" suffixes, which correspond to whether ambiguous
# reads were thrown out.

# Paralog features are collected from a VCF including all known features
# Sequences must be aligned in SAM format to the same reference as was used to generate the VCF

# The total read counts are printed to stdout

# Input files:
#   -Aligned sequences (.SAM)
#   -Variant file (.VCF)

# Output files:
#   -Given a haplotype+suffix the following files are created:
#       .fa    -    a list of the unaligned reads that were estimated to match that haplotype with their bit score
#                   (from -3 to 0) tagged onto their name
#       .txt   -    a feature table for all the aligned reads assigned to the haplotype (with the first row being the
#                   reference haplotype, see far right columns for names/scores)
#   -qualityHist.png: a histogram showing the number of non-missing features per read
#   -exceptionalReads.fa: all low scoring but high "quality" reads are dumped here for troubleshooting purposes. Their
#    paralog matches are listed after the FASTA labels along with score and "quality" (% non-missing features)

# A plot is generated showing the score distributions for each variant, colored by the % of non-missing features in
# each aligned read, Yellow=80-100%, Black=0-20%

from MarkovAligner import *
from collections import defaultdict
from matplotlib import pyplot as plot
from math import floor
import vcf
import pandas as pd
import pysam
import os

# samFileName = "alignedReadsFrankANDYuliaANDMeghan.sam"
samFileName = "margin_align_longer_consensus.sam"
# samFileName = "alignedReadsJasonAdult.sam"
# samFileName = "alignedReadsJasonFetal.sam"
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

featureVectorsByName = dict()
for entry in results:
    featureVectorsByName[entry[0]] = entry[1:]

names = feature_df.index.tolist()[:]        # names of each haplotype (from Alex)
featureArray = feature_df.as_matrix()[:,:]  # set of features for the proposed haplotypes

# ambiguous reads excluded
variantCounts_noTies = {k: 0 for k in names}    # dataset-wide counts of each alignment
variantBins_noTies = {k: list() for k in names} # contains (readName,score) for each aligned read

# ambiguous reads split evenly
variantCounts_all = {k: 0 for k in names}
variantBins_all = {k: list() for k in names}

qualities = list()
for read in results:
    read = list(read[1:])
    # print(len(read))
    quality = 1.0-read.count(None)/float(len(read))

    qualities.append(quality)

model = MarkovAligner()
model.buildModel(featureArray,names,distribution="bernoulliDiscrete")    #generate a model based on the proposed haplotypes


allFeatureVectors = [x[1] for x in results]
nVectors = len(allFeatureVectors)
# nTraining = nVectors*2/3    #unused at the moment

print("Total Feature Vectors: %d"%nVectors)
# print("Training Vectors: %d"%nTraining)

# trainingFeatureVectors = allFeatureVectors[:nTraining]
# testingFeatureVectors = allFeatureVectors[nTraining:]


for i,entry in enumerate(results[:]):
    readName = entry[0]
    readFeatureVector = entry[1:]

    # find the most probable haplotype(s)
    paths,maxPosteriorProb = model.maximumProbabilityPath(readFeatureVector)

    for path in paths:
        haplotypeMatchName = path

        # track alignment counts for each haplotype (splitting ties)
        variantCounts_all[haplotypeMatchName] += 1.0/len(paths)
        variantBins_all[haplotypeMatchName].append((readName,maxPosteriorProb,qualities[i]))

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
                quality = entry[2]

                binScores.append((2**score,quality))

                # featureVector = [x[1:] for x in results if x[0] == entry[0]][0]  #inefficient search
                featureVector = featureVectorsByName[name]
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


featureVectorsSortedByMissingValues = list()
qualityScores = list()

#Create a read dump sorted by a crude measurement of quality (number of missing positions in the feature vector)
for entry in results[:]:
    read = [x if x!=None else 'N' for x in entry[1:]]
    quality = read.count('N')
    qualityScores.append(quality)

    featureVectorsSortedByMissingValues.append((read, quality))

featureVectorsSortedByMissingValues.sort(key=lambda x: x[1])

x = qualityScores

plot.hist(x,25,normed=0,histtype='bar',color=[.118, .565, 1.00],linewidth=0.1)
plot.title("Known Features per Read (quality)")
plot.savefig(outFolder+'/'+"qualityHist.png",dpi=300)

with open(outFolder+'/'+"featureVectorsSortedByMissingValues.txt",'w') as file:
    for entry in featureVectorsSortedByMissingValues:
        file.write('\t'.join(map(str,list(entry[0])))+'\n')

nThresholds = 5
# qualityThresholds = [1.0/nThresholds*n for n in range(nThresholds)]
thresholdedReads = [list() for n in range(nThresholds)]
# print(qualityThresholds)
# colors = [(1.0/nThresholds*n,0.0,(1.0-1.0/nThresholds*n)) for n in range(nThresholds)]
colors = ["#000000","#552B4F","#B93539","#EE4F24","#F99C1E"] #flame LUT

fig,axes = plot.subplots(len(names),1,sharex=True,sharey=True)

# plot the score distributions as log histograms for each haplotype's bin
for h,bin in enumerate(sorted(binnedScoresAll)):
    x = [list() for n in range(nThresholds)]
    for readStats in binnedScoresAll[bin]:
        qualityBinIndex = int(floor(readStats[1]*nThresholds))

        if qualityBinIndex == nThresholds: qualityBinIndex = nThresholds-1  # ensure there is no tailing bin

        x[qualityBinIndex].append(readStats[0])

    for i, thresholdBin in enumerate(x):
        if len(thresholdBin) == 0:
            x[i] = [0]

    # print(len(x))
    # print(len(colors))

    axes[h].hist(x,25,normed=0,histtype='bar',color=colors,linewidth=0.1,range=[0.125,1.0],log=False,stacked=True)
    axes[h].set_ylabel(bin[-3:])

    if h == len(axes)-1:
        axes[h].set_xlabel("Match Likelihood")



plot.show()


