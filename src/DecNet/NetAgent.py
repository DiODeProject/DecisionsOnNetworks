'''
Created on 27 Jan 2017

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

from DecNet.MyTypes import DecisionModel, AgentType, UpdateModel
import math 
import numpy as np
import numpy.random as rand
import scipy.integrate as integrate
import scipy.special
from scipy.stats import norm
import itertools
    
class NetAgent:
    
    DEBUG = False
    
    def __init__(self, agentType, args, debug=False):
        self.opinion = 0
        self.confidence = 0
        self.agentType = agentType
        if (self.agentType == AgentType.SIMPLE):
            self.accuracy = args[0]
        elif (self.agentType == AgentType.DDM):
            self.drift =        args[0]
            self.noiseStdDev =  args[1] 
            self.start =        args[2] # initial position DDM
            self.dt =           args[3]
            self.interrogationT=args[4]
            self.prior =        args[5]
            self.DDMintegration = []
        self.DEBUG = debug
    
    def setPopMeanAndStdDev(self, mean, stdDev):
        self.populationMean = mean
        self.populationStdDev = stdDev
        
    def setBeliefEpsilon(self, eps): self.beliefEpsilon = eps
        
    def setFiniteTimeExponent(self, fte): self.finiteTimeExponent = fte
    
    ## Function computing the perfect-knowledge Log odds 
    def logOddsPerfect(self, decisionVariable):
        pdf = lambda decVar, drift, noise, time: np.exp( (-(decVar - drift * time)**2) / (2 * time * noise * noise) ) / math.sqrt(2 * math.pi * time * noise * noise)
        
        numerator =   pdf( abs(decisionVariable),  abs(self.drift), self.noiseStdDev, self.interrogationT )
        denominator = pdf( abs(decisionVariable), -abs(self.drift), self.noiseStdDev, self.interrogationT )
        #print("num is " +  str(numerator) + " den is " + str(denominator) + " (1-num) = " + str(1-numerator))
        return math.log( numerator/denominator )
    
    ## Function computing the combined Log odds (assuming all agents' drifts are known) 
    def logOddsCombo(self, decisionVariable, all_agents):
        #pdf = lambda decVar, drift, time: np.exp( (-(decVar - drift * time)**2) / (2 * time) ) / math.sqrt(2 * math.pi * time)
        pdf = lambda decVar, drift, noise, time: np.exp( (-(decVar - drift * time)**2) / (2 * time * noise * noise) ) / math.sqrt(2 * math.pi * time * noise * noise)

        numerator = 0
        denominator = 0 
        for agent in all_agents:
            numerator   += pdf( abs(decisionVariable),  abs(agent.drift), self.noiseStdDev, self.interrogationT )
            denominator += pdf( abs(decisionVariable), -abs(agent.drift), self.noiseStdDev, self.interrogationT )
        return math.log( numerator/denominator )

    def integrationFunc(self, acc, decVar, noise, time, sign):
        prob = norm.pdf(acc, self.populationMean, self.populationStdDev)
        if (acc < 0.5): # because agents with accuracy under 0.5 always revert their vote (truncated_low) 
            acc = 1 - acc
        if (sign < 0): # in case of drift with negative signs, it's sufficient to take the mirror accuracy 
            acc = 1 - acc
#         if (sign > 0 and acc < 0.5):
#             acc = 1 - acc
#         if (sign < 0 and acc > 0.5):
#             acc = 1 - acc
        drift = math.sqrt(2) * noise * scipy.special.erfcinv(2 - 2*acc)/math.sqrt(time)
        pdf = np.exp( (-(decVar - drift * time)**2) / (2 * time * noise * noise) ) / math.sqrt(2 * math.pi * time * noise * noise)
        return pdf*prob

    ## Function computing the combined Log odds from distribution (assuming only drift distribution is known) 
    def logOddsDistribution(self, decisionVariable):
        #pdf = lambda drift, decVar, noise, time: np.exp( (-(decVar - drift * time)**2) / (2 * time * noise * noise) ) / math.sqrt(2 * math.pi * time * noise * noise)
        #driftFromAcc = lambda acc, noise, time: math.sqrt(2) * noise * scipy.special.erfcinv(2 - 2*acc)/math.sqrt(time)
        numerator   = integrate.quad( self.integrationFunc, 0, 1, args=( abs(decisionVariable), self.noiseStdDev, self.interrogationT, +1 ) )
        denominator = integrate.quad( self.integrationFunc, 0, 1, args=( abs(decisionVariable), self.noiseStdDev, self.interrogationT, -1 ) )
        #print ("y:" + str(decisionVariable) + " num:" + str(numerator[0]) + " den:" + str(denominator[0]) )
        return math.log( numerator[0]/denominator[0] )
    
    def integrationFuncDriftNormal(self, drift, decVar, noise, time, sign):
        prob = norm.pdf(drift, self.populationMean, self.populationStdDev)
        pdf = np.exp( (-(decVar - (sign*abs(drift)) * time)**2) / (2 * time * noise * noise) ) / math.sqrt(2 * math.pi * time * noise * noise)
        return pdf*prob
    
    def integrationFuncDriftUniform(self, drift, decVar, noise, time, sign):
        pdf = np.exp( (-(decVar - (sign*abs(drift)) * time)**2) / (2 * time * noise * noise) ) / math.sqrt(2 * math.pi * time * noise * noise)
        return pdf

    ## Function computing the combined Log odds from distribution (assuming only drift distribution is known) 
    def estimateConfFromDistribution(self, decisionVariable, time):
#         numerator   = integrate.quad( self.integrationFuncAcc, 0, 1, args=( abs(decisionVariable), self.noiseStdDev, time, +1 ) )
#         denominator = integrate.quad( self.integrationFuncAcc, 0, 1, args=( abs(decisionVariable), self.noiseStdDev, time, -1 ) )
#         print ("y:" + str(decisionVariable) + " num:" + str(numerator[0]) + " den:" + str(denominator[0]) )
#         print("ratio: " + str(numerator[0]/denominator[0]) + " and log: " + str(math.log( numerator[0]/denominator[0] )))
#         numerator   = integrate.quad( self.integrationFuncDriftUniform, self.populationMean-(self.populationStdDev*np.sqrt(12)/2), self.populationMean+(self.populationStdDev*np.sqrt(12)/2), args=( abs(decisionVariable), self.noiseStdDev, time, +1 ) )
#         denominator   = integrate.quad( self.integrationFuncDriftUniform, self.populationMean-(self.populationStdDev*np.sqrt(12)/2), self.populationMean+(self.populationStdDev*np.sqrt(12)/2), args=( abs(decisionVariable), self.noiseStdDev, time, -1 ) )
#         print ("y:" + str(decisionVariable) + " num:" + str(numerator[0]) + " den:" + str(denominator[0]) + " range is [" + str(self.populationMean-(self.populationStdDev*np.sqrt(12)/2)) + "," + str(self.populationMean+(self.populationStdDev*np.sqrt(12)/2)) + "]" )
#         print("ratio: " + str(numerator[0]/denominator[0]) + " and log: " + str(math.log( numerator[0]/denominator[0] )))
        numerator   = integrate.quad( self.integrationFuncDriftNormal, -np.inf, np.inf, args=( abs(decisionVariable), self.noiseStdDev, time, +1 ) )
        denominator   = integrate.quad( self.integrationFuncDriftNormal, -np.inf, np.inf, args=( abs(decisionVariable), self.noiseStdDev, time, -1 ) )
#         print ("y:" + str(decisionVariable) + " num:" + str(numerator[0]) + " den:" + str(denominator[0]) + " range is [" + str(self.populationMean-(self.populationStdDev*np.sqrt(12)/2)) + "," + str(self.populationMean+(self.populationStdDev*np.sqrt(12)/2)) + "]" )
        #print(" num:" + str(numerator[0]) + " den:" + str(denominator[0]) + "ratio: " + str(numerator[0]/denominator[0]) + " and log: " + str(math.log( numerator[0]/denominator[0] )))
        priorComponent = 0 if (self.prior == 0.5) else math.log( self.prior/ (1-self.prior))  # if the sing of decisionVariable is negative, through math.copysign I change the sing of the priorComponent (which is equivalent to power to -1 the log argument)
        return math.log( numerator[0]/denominator[0] ) + math.copysign(priorComponent, decisionVariable)
    
    def logOddsApprox(self, decisionVariable):
        fittedline = np.poly1d([0.1663, 0.5309, 0.1238])
        return( fittedline(abs(decisionVariable)) )
    
    def assignRandomOpinion(self):
        if (rand.random()>0.5):
            return 1
        else:
            return -1
        
    def computeAccuracyFromDrift(self):
        psi = lambda u: np.exp( -u**2 / 2.0 )/ math.sqrt(2*math.pi)
        err = lambda y: integrate.quad( psi, -np.inf, y)
        acc = 1 - err( -self.drift * math.sqrt(self.interrogationT) / self.noiseStdDev )[0]
        return acc
        
    def initialiseOpinion(self, decisionModel, args=None) :
        if (self.agentType == AgentType.SIMPLE):
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
            
        elif (self.agentType == AgentType.DDM):
            ## DDM integration till time T
            t = 0
            self.y = self.start
            while t < self.interrogationT:
                t += self.dt
                self.y += self.drift * self.dt # deterministic step
                self.y += rand.normal(0, self.noiseStdDev) * math.sqrt(self.dt) # stochastic step
                self.DDMintegration.append(self.y)
                
            if (self.y > 0):
                self.opinion = 1
            elif (self.y < 0):
                self.opinion = -1
            
            # compute confidence    
            if (decisionModel == DecisionModel.CONFIDENCE):
                self.accuracy = self.computeAccuracyFromDrift()
                self.confidence = math.log( self.accuracy/(1-self.accuracy) )
            elif (decisionModel == DecisionModel.BELIEF):
                self.accuracy = self.computeAccuracyFromDrift()
                self.confidence = math.log(2.0*self.accuracy)/math.log(2)
            elif (decisionModel == DecisionModel.LOGODDS_PERFECT):
                self.confidence = self.logOddsPerfect(self.y)
            elif (decisionModel == DecisionModel.LOGODDS_COMBO):
                self.confidence = self.logOddsCombo(self.y, args )
            elif (decisionModel == DecisionModel.LOGODDS_DISTRIBUTION):
                #self.confidence = self.logOddsDistribution(self.y)
                self.confidence = self.estimateConfFromDistribution(self.y, self.interrogationT)
            elif (decisionModel == DecisionModel.LOGODDS_APPROX):
                self.confidence = self.logOddsApprox(self.y)
            else:
                self.confidence = 1
                
            if (self.DEBUG): print ("for y=" + str(self.y) + " (and drift:" + str(self.drift) + ") the log-odds conf is " + str(self.confidence))
#             tmp_acc = self.computeAccuracyFromDrift()
#             tmp_conf = math.log( tmp_acc/(1-tmp_acc) )
#             if (self.DEBUG): print ("and from accuracy the confidence would be " + str(tmp_conf))
    
    
    def updateConfidenceOptimQuick(self, all_agents, neighbours, myPreviousOpinion):
        accuracies = []
        opinions = []
        epsilon = 1e-10
        for neigh in neighbours:
            acc = max( epsilon, min( np.exp(all_agents[neigh].confidence) / ( 1 + np.exp(all_agents[neigh].confidence) ), 1-epsilon ))
            accuracies.append( acc ) 
            opinions.append( all_agents[neigh].opinion )
            #print (str(acc) + " -c: " + str(all_agents[neigh].confidence) )
        accuracies.append( self.accuracy )
        opinions.append( myPreviousOpinion )
        #print("accuracies are: " + str(accuracies) )
        #print("opinions are: " + str(opinions) )
        
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
        #print("A1: " + str(self.accuracy) )
        #print("A2: " + str(math.log(comboProbabilityPlus/comboProbabilityNeg)) + " with plus:" + str(comboProbabilityPlus) + " and neg:" + str(comboProbabilityNeg) )
        self.accuracy = max( epsilon, min( 1-epsilon, self.accuracy ))
        if (self.DEBUG): print("Updated accuracy is " + str(self.accuracy))
        if (self.accuracy == 1.0):
            self.confidence = 100
        else:
            self.confidence = math.log( self.accuracy/(1-self.accuracy) )
        #print("C1: " + str(self.confidence) )
        #conf2 = math.log( comboProbabilityPlus/comboProbabilityNeg )
        #acc2 = np.exp(conf2) / ( 1 + np.exp(conf2) )
        #print("acc: " + str(self.accuracy) + " : " + str(acc2) )
        #print("conf:" + str(self.confidence) + " : " + str(conf2) )
#         self.confidence = math.log( comboProbabilityPlus/comboProbabilityNeg )
#         self.accuracy = np.exp(self.confidence) / ( 1 + np.exp(self.confidence) )
#         self.accuracy = max( epsilon, min( 1-epsilon, self.accuracy ))
        
        
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
        
        if (decisionModel == DecisionModel.CONFIDENCE or
            decisionModel == DecisionModel.BELIEF or
            decisionModel == DecisionModel.MAJORITY_RAND or
            decisionModel == DecisionModel.MAJORITY_BIAS or
            decisionModel == DecisionModel.MAJORITY_INHIB or
            decisionModel == DecisionModel.LOGODDS_PERFECT or
            decisionModel == DecisionModel.LOGODDS_COMBO or
            decisionModel == DecisionModel.LOGODDS_DISTRIBUTION or
            decisionModel == DecisionModel.LOGODDS_APPROX):
            ## Weighting my opinion by the confidence
            aggregatedOpinion = self.opinion * self.confidence 
            
            if (updateConfidence == UpdateModel.BELIEF_UP or updateConfidence == UpdateModel.FINITE_TIME):
                oldConf = self.opinion * self.confidence
            
            ## Summing up the weighted opinion of all my neighbours
            if (self.DEBUG):
                if (self.agentType == AgentType.SIMPLE):
                    logAcc = self.accuracy
                elif (self.agentType == AgentType.DDM):
                    logAcc = self.drift
                print('My opinion is ' + str((int)(self.opinion)) + ' (from acc/dirft:' + str(logAcc) + ') and my confidence: ' + 
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
                ## elif (self.method == Method.MAJORITY_BIAS):
                    ## do-nothing
                elif (decisionModel == DecisionModel.MAJORITY_INHIB):
                    self.opinion = 0
            
            ## Updating confidence
            if (updateConfidence == UpdateModel.THETA_UPDATE or updateConfidence == UpdateModel.BELIEF_UP or updateConfidence == UpdateModel.FINITE_TIME):
                self.confidence = min( abs(aggregatedOpinion), 9e+52 )
                if (self.confidence == 0):
                    if (decisionModel == DecisionModel.MAJORITY_RAND or
                        decisionModel == DecisionModel.MAJORITY_BIAS ):
                        self.confidence = 0.5
            elif (updateConfidence == UpdateModel.THETA_NORM):
                self.confidence = abs(aggregatedOpinion)/( len(neighbours)+1 )
            elif (updateConfidence == UpdateModel.OPTIMAL):
                if (decisionModel == DecisionModel.CONFIDENCE):
                    #self.updateConfidenceOptimQuick(all_agents, neighbours, myPreviousOpinion)
                    self.confidence = min( abs(aggregatedOpinion), 9e+52 )
                else:
                    self.confidence = min( abs(aggregatedOpinion), 9e+52 )
                    if (self.confidence == 0):
                        self.confidence = 0.00001
#             elif (updateConfidence == UpdateModel.BELIEF_UP):
#                 oldConf = self.confidence
#                 for neigh in neighbours:
#                     self.confidence += (all_agents[neigh].confidence - oldConf)
            
            if (self.DEBUG and updateConfidence != UpdateModel.NO_UPDATE): print('My new confidence is ' + str(self.confidence))
            
            if (self.DEBUG): print('My new opinion is ' + str(aggregatedOpinion) + ' rounded to ' + str(self.opinion))
        
        elif (decisionModel == DecisionModel.BEST_ACC): ## Warning! The opinion will stop spreading as soon consensus is reached 
            for neigh in neighbours:
                if (self.agentType == AgentType.SIMPLE):
                    if (all_agents[neigh].accuracy > self.accuracy):
                        self.accuracy = all_agents[neigh].accuracy
                        self.opinion = all_agents[neigh].opinion
                elif (self.agentType == AgentType.DDM):
                    if (all_agents[neigh].drift > self.drift) :
                        self.drift = all_agents[neigh].drift
                        self.opinion = all_agents[neigh].opinion

