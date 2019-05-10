'''
Created on 27 Jan 2017

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

from AAMAS2019.MyTypes import DecisionModel, UpdateModel
import math 
import numpy as np
import numpy.random as rand
import itertools
    
class NetAgent:
    
    DEBUG = False
    
    def __init__(self, acc, debug=False):
        self.opinion = 0
        self.confidence = 0
        self.accuracy = acc
        self.DEBUG = debug
    
    def setMeanAccuracyAndStdDev(self, mean, stdDev):
        self.accuracyMean = mean
        self.accuracyStdDev = stdDev
        
    def setBeliefEpsilon(self, eps): self.beliefEpsilon = eps
        
    def setFiniteTimeExponent(self, fte): self.finiteTimeExponent = fte
   
    def assignRandomOpinion(self):
        if (rand.random()>0.5):
            return 1
        else:
            return -1
                
    def initialiseOpinion(self, decisionModel) :
        if (rand.random() < self.accuracy):
            self.opinion = 1
        else:
            self.opinion = -1
        
        # compute confidence
        if (decisionModel == DecisionModel.CONFIDENCE):
            self.confidence = math.log( self.accuracy/(1-self.accuracy) )
            #self.confidence = 1
            #self.accuracy = 0.73105857863
        elif (decisionModel == DecisionModel.BELIEF):
            self.confidence = math.log(2.0*self.accuracy)/math.log(2)
        else:
            self.confidence = 1
    
    
    def updateConfidenceOptimQuick(self, all_agents, neighbours, myPreviousOpinion):
        accuracies = []
        opinions = []
        epsilon = 1e-10
        for neigh in neighbours:
            acc = max( epsilon, min( np.exp(all_agents[neigh].confidence) / ( 1 + np.exp(all_agents[neigh].confidence) ), 1-epsilon ))
            accuracies.append( acc ) 
            opinions.append( all_agents[neigh].opinion )
        accuracies.append( self.accuracy )
        opinions.append( myPreviousOpinion )
        
        ## Computing the probability that the current voting would happen (assuming it's the correct output)
        comboProbabilityPlus = 1
        for n in range(0, len(opinions)):
            if (opinions[n] == self.opinion):
                comboProbabilityPlus *= accuracies[n]
            else:
                comboProbabilityPlus *= (1 - accuracies[n])
        
        ## Computing the probability that the opposite voting would happen (assuming it's the incorrect output) 
        comboProbabilityNeg = 1
        for n in range(0, len(opinions)):
            if (opinions[n] == self.opinion):
                comboProbabilityNeg *= (1 - accuracies[n])
            else:
                comboProbabilityNeg *= accuracies[n]
            
        ## Computing how likely is the opposite aggregate opinion is correct
        self.accuracy = comboProbabilityPlus/(comboProbabilityPlus + comboProbabilityNeg)
        self.accuracy = max( epsilon, min( 1-epsilon, self.accuracy ))
        if (self.DEBUG): print("Updated accuracy is " + str(self.accuracy))
        if (self.accuracy == 1.0):
            self.confidence = 100
        else:
            self.confidence = math.log( self.accuracy/(1-self.accuracy) )

        
        
    def updateConfidenceOptimFull(self, all_agents, neighbours):
        # invert confidence-weighting function to calculate expected neighbour accuracy
        accuracies = []
        confidences = []
        for neigh in neighbours:
            accuracies.append( np.exp(all_agents[neigh].confidence) / ( 1 + np.exp(all_agents[neigh].confidence) ) )
            confidences.append( all_agents[neigh].confidence )
        accuracies.append( self.accuracy )
        confidences.append( self.confidence )
        #print("accuracies are: " + str(accuracies) )
        
        updatedAcc = 0
        allCombinations = list( itertools.product( [-1,1], repeat=len(accuracies)) )  # all vote combinations
        if (self.DEBUG): print("Running combinatorics update with " + str(len(allCombinations)) + " combinations")
        for combo in allCombinations:
            weightedVote = 0
            for n in range(0, len(combo)): # combining votes with given confidences 
                weightedVote += combo[n] * confidences[n]
            #print("For combo " + str(combo) + " the weightedVote is " + str(weightedVote) )
            if weightedVote > 0: # check if the combined votes are Positive (we increment accuracies assuming positive vote, it would give the same results by using negative)
                comboProbability = 1
                for n in range(0, len(combo)):
                    if (combo[n] == 1):
                        comboProbability *= accuracies[n]
                    else:
                        comboProbability *= (1 - accuracies[n])
                updatedAcc += comboProbability
                #print("And comboProb is " + str(comboProbability) + " and upAcc: " + str(updatedAcc))
        
        self.accuracy = updatedAcc
        if (self.DEBUG): print("Updated accuracy is " + str(updatedAcc))
        self.confidence = math.log( self.accuracy/(1-self.accuracy+0.00000000001) )
        
    
    def integrateInfoFromNeighbours(self, decisionModel, all_agents, neighbours, updateConfidence):
        ## Weighting my opinion by the confidence
        aggregatedOpinion = self.opinion * self.confidence 
        
        if (updateConfidence == UpdateModel.BELIEF_UP or updateConfidence == UpdateModel.FINITE_TIME):
            oldConf = self.opinion * self.confidence
        
        ## Summing up the weighted opinion of all my neighbours
        if (self.DEBUG):
            print('My opinion is ' + str((int)(self.opinion)) + ' (from accuracy :' + str(self.accuracy) + ') and my confidence: ' + 
                          str(self.confidence) + ' combined to: ' + str(aggregatedOpinion))
        for neigh in neighbours:
            if (self.DEBUG): print('Receive from neigh ' + str(neigh) + ' its opinion (' +
                              str((int)(all_agents[neigh].opinion)) + ') and its confidence: ' + str(all_agents[neigh].confidence) +
                              ' (' + str(all_agents[neigh].opinion * all_agents[neigh].confidence) + ')')
            if (updateConfidence == UpdateModel.BELIEF_UP):
                belief_value = 1/(1+len(neighbours)) if (self.beliefEpsilon == 0) else self.beliefEpsilon
                aggregatedOpinion += ( ( (all_agents[neigh].opinion*all_agents[neigh].confidence) - oldConf) * belief_value )  #/(len(all_agents)-1) )
            elif (updateConfidence == UpdateModel.FINITE_TIME):
                tmp_diff = (all_agents[neigh].opinion*all_agents[neigh].confidence) - oldConf
                belief_value = 1/(1+len(neighbours)) if (self.beliefEpsilon == 0) else self.beliefEpsilon
                #aggregatedOpinion += self.beliefEpsilon * np.sign(tmp_diff) * (np.abs(tmp_diff)**self.beliefEpsilon)
                aggregatedOpinion += belief_value * np.sign(tmp_diff) * (np.abs(tmp_diff)**self.finiteTimeExponent)
                if abs(aggregatedOpinion) > 9e+52: break              
            else:
                aggregatedOpinion += (all_agents[neigh].opinion * all_agents[neigh].confidence)
        
    
        ## Determining the new opinion and the new confidence
        #myPreviousOpinion = self.opinion 
        if (aggregatedOpinion>0):
            self.opinion = 1
        elif (aggregatedOpinion<0):
            self.opinion = -1
        elif (aggregatedOpinion==0):
            if (decisionModel == DecisionModel.CONFIDENCE or 
                decisionModel == DecisionModel.BELIEF):
                self.opinion = self.assignRandomOpinion()
            elif (decisionModel == DecisionModel.MAJORITY_RAND):
                self.opinion = self.assignRandomOpinion()
        
        ## Updating confidence
        if (updateConfidence == UpdateModel.BELIEF_UP or updateConfidence == UpdateModel.FINITE_TIME):
            self.confidence = min( abs(aggregatedOpinion), 9e+52 )
            if (self.confidence == 0):
                if (decisionModel == DecisionModel.MAJORITY_RAND):
                    self.confidence = 0.5
        elif (updateConfidence == UpdateModel.OPTIMAL):
            if (decisionModel == DecisionModel.CONFIDENCE):
                #self.updateConfidenceOptimQuick(all_agents, neighbours, myPreviousOpinion)
                self.confidence = min( abs(aggregatedOpinion), 9e+52 )
            else:
                self.confidence = min( abs(aggregatedOpinion), 9e+52 )
                if (self.confidence == 0):
                    self.confidence = 0.00001
        
        if (self.DEBUG and updateConfidence != UpdateModel.NO_UPDATE): print('My new confidence is ' + str(self.confidence))        
        if (self.DEBUG): print('My new opinion is ' + str(aggregatedOpinion) + ' rounded to ' + str(self.opinion))
        

