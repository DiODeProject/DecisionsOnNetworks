'''
Created on 15 Mar 2019

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

from AsynchKicks.MyTypes import AgentType, UpdateModel
import math 
import numpy as np
import numpy.random as rand
import scipy.integrate as integrate
import scipy.special
from scipy.stats import norm
import itertools
    
class DdmAgent:
    
    DEBUG = False
    
    def __init__(self, myid, agentType, updateModel, args, debug=False):
        self.myid = myid
        self.opinion = 0
        self.confidence = 0
        self.agentType = agentType
        self.updateModel = updateModel
        if (self.agentType == AgentType.SIMPLE):
            self.accuracy = args[0]
        elif (self.agentType == AgentType.DDM):
            self.drift =        args[0]
            self.noiseStdDev =  args[1] 
            self.start =        args[2] # initial position DDM
            self.dt =           args[3]
            self.threshold =    args[4]
            self.prior =        args[5]
            self.DDMintegration = []
            self.DDMintegration_time = []
        self.DEBUG = debug
    
    def setMeanAndStdDev(self, mean, stdDev):
        self.popMean = mean
        self.popStdDev = stdDev
        
    def setBeliefEpsilon(self, eps): self.beliefEpsilon = eps
        
    def setFiniteTimeExponent(self, fte): self.finiteTimeExponent = fte
    
    ## Function computing the perfect-knowledge Log odds 
    def logOddsPerfect(self, decisionVariable, time):
        pdf = lambda decVar, drift, noise, time: np.exp( (-(decVar - drift * time)**2) / (2 * time * noise * noise) ) / math.sqrt(2 * math.pi * time * noise * noise)
        
        numerator =   pdf( abs(decisionVariable),  abs(self.drift), self.noiseStdDev, time )
        denominator = pdf( abs(decisionVariable), -abs(self.drift), self.noiseStdDev, time )
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

    def integrationFuncAcc(self, acc, decVar, noise, time, sign):
        prob = norm.pdf(acc, self.popMean, self.popStdDev)
        if (acc < 0.5): # because agents with accuracy under 0.5 always revert their vote (truncated_low) 
            acc = 1 - acc
        if (sign < 0): # in case of drift with negative signs, it's sufficient to take the mirror accuracy 
            acc = 1 - acc
        drift = math.sqrt(2) * noise * scipy.special.erfcinv(2 - 2*acc)/math.sqrt(time)
        pdf = np.exp( (-(decVar - drift * time)**2) / (2 * time * noise * noise) ) / math.sqrt(2 * math.pi * time * noise * noise)
        return pdf*prob

    ## Function computing the combined Log odds from distribution (assuming only drift distribution is known) 
    def logOddsDistribution(self, decisionVariable):
        #pdf = lambda drift, decVar, noise, time: np.exp( (-(decVar - drift * time)**2) / (2 * time * noise * noise) ) / math.sqrt(2 * math.pi * time * noise * noise)
        #driftFromAcc = lambda acc, noise, time: math.sqrt(2) * noise * scipy.special.erfcinv(2 - 2*acc)/math.sqrt(time)
        numerator   = integrate.quad( self.integrationFuncAcc, 0, 1, args=( abs(decisionVariable), self.noiseStdDev, self.interrogationT, +1 ) )
        denominator = integrate.quad( self.integrationFuncAcc, 0, 1, args=( abs(decisionVariable), self.noiseStdDev, self.interrogationT, -1 ) )
        #print ("y:" + str(decisionVariable) + " num:" + str(numerator[0]) + " den:" + str(denominator[0]) )
        return math.log( numerator[0]/denominator[0] )

    def integrationFuncDriftNormal(self, drift, decVar, noise, time, sign):
        prob = norm.pdf(drift, self.popMean, self.popStdDev)
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
#         numerator   = integrate.quad( self.integrationFuncDriftUniform, self.popMean-(self.popStdDev*np.sqrt(12)/2), self.popMean+(self.popStdDev*np.sqrt(12)/2), args=( abs(decisionVariable), self.noiseStdDev, time, +1 ) )
#         denominator   = integrate.quad( self.integrationFuncDriftUniform, self.popMean-(self.popStdDev*np.sqrt(12)/2), self.popMean+(self.popStdDev*np.sqrt(12)/2), args=( abs(decisionVariable), self.noiseStdDev, time, -1 ) )
#         print ("y:" + str(decisionVariable) + " num:" + str(numerator[0]) + " den:" + str(denominator[0]) + " range is [" + str(self.popMean-(self.popStdDev*np.sqrt(12)/2)) + "," + str(self.popMean+(self.popStdDev*np.sqrt(12)/2)) + "]" )
#         print("ratio: " + str(numerator[0]/denominator[0]) + " and log: " + str(math.log( numerator[0]/denominator[0] )))
        numerator   = integrate.quad( self.integrationFuncDriftNormal, -np.inf, np.inf, args=( abs(decisionVariable), self.noiseStdDev, time, +1 ) )
        denominator   = integrate.quad( self.integrationFuncDriftNormal, -np.inf, np.inf, args=( abs(decisionVariable), self.noiseStdDev, time, -1 ) )
#         print ("y:" + str(decisionVariable) + " num:" + str(numerator[0]) + " den:" + str(denominator[0]) + " range is [" + str(self.popMean-(self.popStdDev*np.sqrt(12)/2)) + "," + str(self.popMean+(self.popStdDev*np.sqrt(12)/2)) + "]" )
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
        
    def computeAccuracyFromDrift(self, time):
        psi = lambda u: np.exp( -u**2 / 2.0 )/ math.sqrt(2*math.pi)
        err = lambda y: integrate.quad( psi, -np.inf, y)
        acc = 1 - err( -self.drift * math.sqrt(time) / self.noiseStdDev )[0]
        return acc
    
    def initialiseDDM(self):
        self.y = self.start # initialise the DDM
        self.DDMintegration.append(self.y)
        self.DDMintegration_time.append(0)
    
    def updateOpinion(self):
        if self.y >= self.threshold:
            self.opinion = +1
        elif self.y <= -self.threshold:
            self.opinion = -1
    
    def integrateIndividualEvidence(self, _time) :
        self.y += self.drift * self.dt # deterministic step
        self.y += rand.normal(0, self.noiseStdDev) * math.sqrt(self.dt) # stochastic step
        self.DDMintegration.append(self.y)
        self.DDMintegration_time.append(_time)
        #print("y val: " + str(self.y))
        self.updateOpinion()
        #if self.opinion != 0: print("My perfect conf (a:" + str(self.myid) + ") at time " + str(_time) + " would be " + str(self.logOddsPerfect(self.y, _time)))
#         if self.opinion != 0: 
#             acc = self.computeAccuracyFromDrift(_time)
#             print("just for test: " + str( math.log( acc/(1-acc)) ) )
        return self.opinion

    def kick(self, kickerOpinion, time) :
        thresh = kickerOpinion * self.threshold
        if (self.updateModel == UpdateModel.NO_UPDATE):
            kicksize = 0
            return self.opinion
        if (self.updateModel == UpdateModel.THRESH_KICK):
            kicksize = thresh
        if (self.updateModel == UpdateModel.CONF_KICK):
            priorComponent = 0 if (self.prior == 0.5) else math.log( self.prior/ (1-self.prior))  # if the sing of kickerOpinion is negative, through math.copysign I change the sing of the priorComponent (which is equivalent to power to -1 the log argument)
            kicksize = kickerOpinion * self.estimateConfFromDistribution(thresh, time) - math.copysign(priorComponent, kickerOpinion) 
        self.y += kicksize 
        #print("kicksize is " + str(kicksize) + " which put my y to " + str(self.y))
        self.DDMintegration.append(self.y)
        self.DDMintegration_time.append(time)
        self.updateOpinion()
        return self.opinion
    
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
        
