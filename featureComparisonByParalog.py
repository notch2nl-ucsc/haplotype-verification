#!/usr/bin/env python2.7

# This script takes a VCF and identifies unique features for each NOTCH2 or NOTCH2NL paralog assuming the name contains
# at some position "NOTCH2NL-P" or "NOTCH-P" where P is the paralog.

# Input:
#   VCF folder and filename (use empty string if VCF in same folder)
#   Suffix to tag the output with. Base output name is the VCF filename

# Output:
#   A text file containing a table of unique features and their ONE-BASED positions in the consensus used for VCF


import vcf
from collections import defaultdict
import pandas as pd
import copy


folder = "UpdatedVariants_6-1-17"
# vcfFilename = "all_tx_variants_minusPos0.vcf"
vcfFilename = "transcripts.vcf"
# vcfFilename = "variants_longer.vcf"
outputFileSuffix = "_test"  #change this if you don't want to overwrite the last run

#----------------------------------------------------------------------------------------------------------------------#
# from IAN FIDDES

recs = list(vcf.Reader(open(folder+'/'+vcfFilename)))
# recs = [x for x in recs if x.POS <= 1100]

intervals = []

for r in recs:
    start, stop = r._compute_coordinates_for_indel()
    intervals.append([start - 1, stop, r])  # convert to 0-based

# how many distinct feature vectors are there?
feature_vectors = defaultdict(list)  # why defaultdict? <- intializes all values as an empty list
for r in recs:
    for sample in r.samples[:-1]:  # exclude consensus
        gt = sample.gt_bases.split('/')[0]  # remove duplicate SNP (homozygous)
        feature_vectors[sample.sample].append(gt)  # store sample name:[SNP1,SNP2...SNPn] in dictionary

feature_df = pd.DataFrame.from_dict(feature_vectors).T
# feature_df = feature_df.drop_duplicates()

#----------------------------------------------------------------------------------------------------------------------#
variantStartPositions_oneBased = [interval[0]+1 for interval in intervals]

fileSuffix = vcfFilename.split('.')[0]

names = feature_df.index.tolist()[:]        # names of each haplotype (from Alex)
featureArray = feature_df.as_matrix()[:,:]  # set of features for the proposed haplotypes


def extractUniqueFeatures(names,featureArray,paralogSubset=None,outputAsList=True):
    '''
    From a list of features from NOTCH variants, produce a set of unique features at each position across paralogs.
    <names> is the list of variants, and <featureArray> is their extracted SNPs.

    To compare a subset of paralogs, use paralogSubset=set(['A','B'])

    Output is a list of vectors for each paralog, containing a list object for each variant position. Empty
    lists indicate no unique features at a given position.

    outputAsList==False will retain internal sets in the array at each position, ideal for further comparison...
    '''

    paralogFeatures = dict()
    for name,vector in zip(names,featureArray):
        # determine the NOTCH paralog type from VCF file
        if "NOTCH2NL-" in name.upper():
            name = name.split('-')
            paralog = None
            for i,fragment in enumerate(name):
                if fragment == "NOTCH2NL":
                    paralog = name[i+1]

            # print("%s=%s"%(name,paralog))
            valid = True
        elif "NOTCH2-" in name.upper():
            paralog = 'N2'
            # print("%s=%s"%(name,paralog))
            valid = True
        else:
            print("WARNING: non NOTCH2 paralog in VCF: %s"%name)
            valid = False

        if valid:
            if paralogSubset != None:
                if paralog in paralogSubset:
                    # initialize paralog feature vector
                    if paralog not in paralogFeatures:
                        paralogFeatures[paralog] = [set() for feature in vector]

                    # add features if they are of the same paralog
                    for i,feature in enumerate(paralogFeatures[paralog]):
                        feature.add(vector[i])

            else:
                # initialize paralog feature vector
                if paralog not in paralogFeatures:
                    paralogFeatures[paralog] = [set() for feature in vector]

                # add features if they are of the same paralog
                for i, feature in enumerate(paralogFeatures[paralog]):
                    feature.add(vector[i])


    paralogUniqueFeatures = dict()

    for paralog1 in paralogFeatures:
        queryParalogFeatures = copy.deepcopy(paralogFeatures[paralog1])

        for paralog2 in paralogFeatures:
            if paralog2 != paralog1:
                for i,feature in enumerate(paralogFeatures[paralog2]):
                    queryParalogFeatures[i] = queryParalogFeatures[i] - feature

        if outputAsList:
            paralogUniqueFeatures[paralog1] = list(map(list,queryParalogFeatures))
        else:
            paralogUniqueFeatures[paralog1] = queryParalogFeatures

    return paralogUniqueFeatures


paralogUniqueFeatures = extractUniqueFeatures(names,featureArray,set(['A','B']))

with open("paralogUniqueFeatures"+outputFileSuffix+".txt",'w') as file:
    file.write('\t'+'\t'.join(list(map(str, variantStartPositions_oneBased)))+'\n')

    for paralog in sorted(paralogUniqueFeatures):
        featureVectorPrintable = list()

        for feature in paralogUniqueFeatures[paralog]:
            if len(feature) > 1:
                feature = [','.join(feature)]

            if len(feature) > 0:
                featureVectorPrintable.append(feature[0])
            else:
                featureVectorPrintable.append('-')
        file.write('\t'.join([paralog]+featureVectorPrintable)+'\n')


