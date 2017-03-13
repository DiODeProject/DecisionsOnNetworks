'''
Created on 27 Jan 2017

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

from DecNet.MyTypes import DecisionModel 
from DecNet.MyTypes import AgentType 
import math 
import numpy as np
import numpy.random as rand
import scipy.integrate as integrate
    
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
            self.DDMintegration = []
        self.DEBUG = debug
            
    ## Function computing the Log odds 
    def logOdds(self, decisionVariable, all_agents):
        pdf = lambda decVar, drift, time: np.exp( (-(decVar - drift * time)**2) / (2 * time) ) / math.sqrt(2 * math.pi * time)
#         numerator = pdf(decisionVariable, self.drift, self.interrogationT)
        numerator = 0
        denominator = 0
        for agent in all_agents:
#             numerator   += pdf(abs(decisionVariable),  abs(agent.drift), self.interrogationT)
#             denominator += pdf(abs(decisionVariable), -abs(agent.drift), self.interrogationT)
            numerator   += pdf(decisionVariable,  agent.drift, self.interrogationT)
            denominator += pdf(decisionVariable, -agent.drift, self.interrogationT)
        return numerator/denominator
    
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
            else:
                self.confidence = 1
            
        elif (self.agentType == AgentType.DDM):
            ## DDM integration till time T
            t = 0
            self.y = self.start
            while t <= self.interrogationT:
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
            elif (decisionModel == DecisionModel.LOGODDS):
                self.confidence = self.logOdds(self.y, args)
            else:
                self.confidence = 1
            if (self.DEBUG): print ("for y=" + str(self.y) + " (and drift:" + str(self.drift) + ") the log-odds conf is " + str(self.confidence))
    
    def integrateInfoFromNeighbours(self, decisionModel, all_agents, neighbours):
        
        if (decisionModel == DecisionModel.CONFIDENCE or
            decisionModel == DecisionModel.MAJORITY_RAND or
            decisionModel == DecisionModel.MAJORITY_BIAS or
            decisionModel == DecisionModel.MAJORITY_INHIB or
            decisionModel == DecisionModel.LOGODDS):
            ## Weighting my opinion by the confidence
            aggregatedOpinion = self.opinion * self.confidence 
            
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
                aggregatedOpinion += (all_agents[neigh].opinion * all_agents[neigh].confidence)
        
            ## Determining the new opinion and the new confidence
            if (aggregatedOpinion>0):
                self.opinion = 1
            elif (aggregatedOpinion<0):
                self.opinion = -1
            elif (aggregatedOpinion==0):
                if (decisionModel == DecisionModel.CONFIDENCE):
                    self.opinion = self.assignRandomOpinion()
                elif (decisionModel == DecisionModel.MAJORITY_RAND):
                    self.opinion = self.assignRandomOpinion()
                ## elif (self.method == Method.MAJORITY_BIAS):
                    ## do-nothing
                elif (decisionModel == DecisionModel.MAJORITY_INHIB):
                    self.opinion = 0
                    
            if (decisionModel == DecisionModel.CONFIDENCE):
                self.confidence = abs(aggregatedOpinion)
            
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

