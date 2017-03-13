'''
Created on 10 Oct 2016

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

# import math
import numpy.random as rand
import math 
from enum import Enum

class Method(Enum):
    CONFIDENCE = 1
    MAJORITY_RAND = 2
    MAJORITY_BIAS = 3
    MAJORITY_INHIB= 4
    CDCI = 5
        
class DDM:
    converged = 0
    positive = 0
    negative = 0
    
    def __init__(self, drift, stdDev, start, dt, thresh):
        self.drift = drift
        self.stdDev = stdDev
        self.y = start
        self.dt = dt
        self.decided = False
        self.confidence = 0
        self.threshold = thresh
        self.inbox = []
        self.method = Method.CONFIDENCE
        
    def reset(self, start):
        self.y = start
        self.decided = False
        self.confidence = 0
        self.inbox = []
    
    @staticmethod
    def resetGlobalCounters(ddms):
        DDM.positive = 0
        DDM.negative = 0
        for ddm in ddms:
            if (ddm.y >= ddm.threshold):
                DDM.positive += 1
            if (ddm.y <= -ddm.threshold):
                DDM.negative += 1
    
    def step(self):
        if (not self.decided) :
            self.y += self.drift * self.dt # deterministic step
            self.y += rand.normal(0, self.stdDev) * math.sqrt(self.dt) # stochastic step
            if (self.y >= self.threshold or self.y <= -self.threshold):
                self.decided = True
                DDM.converged += 1
                if (self.y >= self.threshold):
                    self.y = 1
                    DDM.positive += 1
                if (self.y <= -self.threshold):
                    self.y = -1
                    DDM.negative += 1
                return True
        return False
    
    ## Function computing the Log odds correct from the convergence time
    def logOdds(self, convTime):
        #return 1/convTime
        
        ## Linear fit done over the function gright from the Farkas' paper (eq.8) using:
        ## 20 drifts in [0,1], L=2, alpha=0.5, sigma=1, T=20
        ## ArXiv's formula: Conf = 0.765574 - 0.0139741 t
        ## Journal formula: Conf = 0.494736 - 0.0125872 t
        ## Thomas's corrected version: Conf = 0.108083 - 0.00127828 t
        return (0.108083 - 0.00127828 * convTime)
        ## 300 drifts in N[0.1,0.3], L=2, alpha=0.5, sigma=0.5, T=10
        #return (0.0734047 - 0.000747686 * convTime)
    
    def assignRandomOpinion(self):
        if (rand.random()>0.5):
            return 1
        else:
            return -1
    
    def integrateInfoFromNeighbours(self, tmp_ddms, neighbours, DEBUG):
        ## Temporarily removing the agent from the converged count for then reintroducing it at the end of the update
        if (self.y > 0):
            DDM.positive -= 1
        elif (self.y < 0):
            DDM.negative -= 1

        if (self.method == Method.CDCI):
            neighList = list(neighbours)
            rand.shuffle(neighList)
            for neigh in neighList:
                if (tmp_ddms[neigh].y == 0):
                    continue
                # compute the update probability from the difference between mine and his confidence 
                confDiff = tmp_ddms[neigh].confidence - self.confidence
                maxConf = self.logOdds(0) # store the max-confidence value 
                confDiff = (confDiff + maxConf)/(2*maxConf) # normalize confDiff
                #print( str(confDiff) + " from " + str(self.confidence) + " and " + str(tmp_ddms[neigh].confidence))
                # perform the update  
                if (rand.random() < confDiff):
                    if (self.y == 0): ## RECRUITMENT
                        self.y = tmp_ddms[neigh].y
                        self.confidence = tmp_ddms[neigh].confidence
                    elif (self.y != 0 and self.y != tmp_ddms[neigh].y): ## CROSS-INHIBITION
                        self.y = 0
                        self.confidence = 0
#                     elif (self.y != 0 and self.y == tmp_ddms[neigh].y): ## Just updating my confidence
#                         print(" I was " + str(self.confidence) + " ad neigh " + str(tmp_ddms[neigh].confidence))
#                         self.confidence = min(maxConf, (self.confidence + tmp_ddms[neigh].confidence*0.005))
#                         print(" I am " + str(self.confidence) )
                break
        else:
            ## Weighting my opinion by the logOdds metric
            myOpinion = self.y * self.confidence 
            #self.logOdds(self.convTime)
            
            ## Summing up the weighted opinion of all my neighbours
            neighsOpinion = 0
            if (DEBUG): print('My opinion is ' + str((int)(self.y)) + ' and my confidence: ' + 
                              str(self.confidence) + ' combined to: ' + str(myOpinion))
            for neigh in neighbours:
                if (DEBUG): print('Receive from neigh ' + str(neigh) + ' its opinion (' +
                                  str((int)(tmp_ddms[neigh].y)) + ') and its confidence: ' + str(tmp_ddms[neigh].confidence) +
                                  ' (' + str(tmp_ddms[neigh].y * tmp_ddms[neigh].confidence) + ')')
                neighsOpinion += (tmp_ddms[neigh].y * tmp_ddms[neigh].confidence)
            if (DEBUG): print('Others opinion is ' + str(neighsOpinion))
            
            ## Integrating my opinion with the neighbours's opinion
            myOpinion = myOpinion + neighsOpinion
        
            ## Determining the new opinion and the new confidence
            if (myOpinion>0):
                self.y = 1
            elif (myOpinion<0):
                self.y = -1
            elif (myOpinion==0):
                if (self.method == Method.CONFIDENCE):
                    self.y = self.assignRandomOpinion()
                elif (self.method == Method.MAJORITY_RAND):
                    self.y = self.assignRandomOpinion()
                ## elif (self.method == Method.MAJORITY_BIAS):
                    ## do-nothing
                elif (self.method == Method.MAJORITY_INHIB):
                    self.y = 0
            
            if (self.method == Method.CONFIDENCE):
                self.confidence = abs(myOpinion)
            
            if (DEBUG): print('My new opinion is ' + str(myOpinion) + ' rounded to ' + str(self.y))
                    
        if (self.y > 0):
            DDM.positive += 1
        elif (self.y < 0):
            DDM.negative += 1  
            
    