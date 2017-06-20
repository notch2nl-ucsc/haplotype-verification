#!/usr/bin/env python2.7

# Markov Aligner by Ryan Lorig-Roach (rlorigro@ucsc.edu)

import random
import numpy
import sys

class MarkovAligner:
    '''
    An extremely simple HMM designed for sequences that are already aligned to a consensus. No sequence length variation is
    possible, and no recombination is allowed between paths, which makes probability calculations trivial.
    '''
    def __init__(self,bernoulliProbability=0.9):
        self.paths = dict()
        self.path = list()
        self.sequence = None
        self.bernoulliProbability = bernoulliProbability
        self.distribution = "bernoulliDiscrete"
        self.distributions = {"bernoulliDiscrete":self.bernoulliDiscreteDistribution,
                              "mixedBernoulliDiscrete":self.mixedBernoulliDiscreteDistribution}

    class mixedBernoulliDiscreteDistribution:
        '''
        A distribution which has multiple possible match values: match or non-match with probability p or q=(1-p).
        Each node in the model contains this distribution, which a match or "symbol" parameter specific to a nucleotide
        or feature at a given position from a template sequence.
        '''

        def __init__(self,symbol=set(),p=0.0,antisymbol='*'):
            self.symbol = symbol    #to be matched (typically a set of sequences or characters)
            self.antisymbol = antisymbol    #all non matches for this state are returned as this during sampling
            self.p = p  #probability of the symbol given this state

            if len(self.symbol)==0:
                self.p=0.9

        def probability(self,symbol):

            if symbol in self.symbol:
                return numpy.log2(self.p)
            else:
                # if symbol!= None:
                #     print('\n')
                #     print(symbol)
                #     print(self.symbol)
                return numpy.log2(1-self.p)

        def sample(self):   #currently unused (no sampling method exists in MarkovAligner class)

            randomFloat = random.random()

            if randomFloat < self.p:
                return self.symbol.pop()
            else:
                return self.antisymbol

    class bernoulliDiscreteDistribution:
        '''
        A distribution which has only two possible values: match or non-match with probability p or q=(1-p).
        Each node in the model contains this distribution, which a match or "symbol" parameter specific to a nucleotide
        or feature at a given position from a template sequence.
        '''

        def __init__(self,symbol,p,antisymbol='*'):
            self.symbol = symbol    #to be matched (typically a sequence or character)
            self.antisymbol = antisymbol    #all non matches for this state are returned as this during sampling
            self.p = p  #probability of the symbol given this state

        def probability(self,symbol):

            if symbol == self.symbol:
                return numpy.log2(self.p)
            else:
                return numpy.log2(1-self.p)

        def sample(self):   #currently unused (no sampling method exists in MarkovAligner class)

            randomFloat = random.random()

            if randomFloat < self.p:
                return self.symbol
            else:
                return self.antisymbol

    def buildModel(self,sequences,pathNames=None,probabilities=None,distribution="bernoulliDiscrete"):
        '''
        Build a model using a set of template sequences, producing one independent path for each template
        '''

        self.distribution = distribution

        if pathNames == None:
            pathNames = list(map(str,range(0,len(sequences))))

        for s,sequence in enumerate(sequences):
            name = pathNames[s]

            self.paths[name] = list()
            self.path = self.paths[name]

            self.buildPath(sequence,probabilities=probabilities)

        # self.clear()

    def buildPath(self,sequence,probabilities=None,):
        '''
        Helper function for building the model
        '''
        if probabilities == None:
            probabilities = [0.9]*len(sequence)

        # print(sequence)

        for s,symbol in enumerate(sequence):
            probability = probabilities[s]
            # print(symbol)

            self.path.append(self.distributions[self.distribution](symbol, probability))

    def clear(self):
        self.path = None
        self.sequence = None

    def forward(self,sequence,pathName):
        '''
        "Forward algorithm" for a linear HMM with no bifurcation. Find the product of probabilities P(Xi|Pi) for each
        observed feature Xi at its i-th position in the path.
        '''

        a = zip(self.paths[pathName],sequence)

        probabilities = [x.probability(y) for x,y in a]

        return sum(probabilities)   #sum of logs is product of probabilities

    def maximumProbabilityPath(self,sequence):
        '''
        Find the maximum probability path (given the full model) for an observed sequence, and return possible ties with
        the score P(maxPath)/sum(P(allPaths))
        '''
        results = list()
        pSum = None

        for pathName in self.paths:
            p = self.forward(sequence,pathName)

            results.append((pathName,p))

            if pSum != None:
                pSum = numpy.logaddexp2(p,pSum)     # use log2 addition formula to prevent underflow
            else:
                pSum = p

        results.sort(key=lambda x: x[1], reverse=True)
        viterbiPaths = [x[0] for x in results if x[1] == results[0][1]]     # find all top-scoring paths

        pMax = results[0][1]    # the top score in sorted results

        maximumProbability = pMax - pSum  # log probability (bits)

        return viterbiPaths,maximumProbability



# # Testing
# sequences = [['a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a']*900,
#              ['b','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b','b']*900,
#              ['c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c']*900]
#
# names = ['A','B','C']
#
# model = MarkovAligner()
# model.buildModel(sequences,names)
#
# testSeq = ['d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d','d']*900
#
# p = model.forward(testSeq,'A')
# q = model.forward(testSeq,'B')
#
# print(2**p)
# print(2**q)
#
# print(model.maximumProbabilityPath(testSeq))

# # More testing
# sequences = [[set(['a'])]*40,
#              [set(['a','b'])]*40,
#              [set(['c'])]*40]
#
# names = ['A','B','C']
#
# model = MarkovAligner()
# model.buildModel(sequences,names,distribution="mixedBernoulliDiscrete")
#
# testSeq = ['d']*40
#
# p = model.forward(testSeq,'A')
# q = model.forward(testSeq,'B')
#
# print(2**p)
# print(2**q)
#
# print(model.maximumProbabilityPath(testSeq))
