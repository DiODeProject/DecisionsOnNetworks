'''
Created on 09 Mar 2017

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

import scipy.special
import math
import numpy as np
import numpy.random as rand
from DecNet.NetAgent import NetAgent
from DecNet.MyTypes import AgentType 
from DecNet.MyTypes import DecisionModel 

def computeDriftFromAccuracy(accuracy, noise, interrogationTime):
        return math.sqrt(2) * noise * scipy.special.erfcinv(2 - 2*accuracy)/math.sqrt(interrogationTime)

if __name__ == '__main__':
    accuracy = 0.8217
    noiseStdDev = 0.5
    interrogationTime = 3
    DDMstart = 0
    dt = 0.01
    seed = 221189
    runs = 10000
    
    rand.seed(seed)
    
    #     acc = rand.normal(self.accuracyMean, self.accuracyStdDev)
#     while (acc < 0 or acc >=  1):
#         acc = rand.normal(self.accuracyMean, self.accuracyStdDev)
    driftRate = computeDriftFromAccuracy(accuracy, noiseStdDev, interrogationTime)
    print("Drift rate is: " + str(driftRate))
    
    args = []
    args.append( driftRate ) 
    args.append( noiseStdDev )
    args.append( DDMstart )
    args.append( dt )
    args.append( interrogationTime )
     
    agent = NetAgent(AgentType.DDM, args, True)
    acc = agent.computeAccuracyFromDrift()
    print("Accuracy is: " + str(acc))
    
    agents = []
    for i in range(0,10):
        args[0] = args[0] + i/100
        agents.append( NetAgent(AgentType.DDM, args, False) )
    
    pdf = lambda decVar, drift, time: np.exp( (-(decVar - drift * time)**2) / (2 * time) ) / math.sqrt(2 * math.pi * time)
    print( pdf(0.4, 0.1, 3) )
    
    success = 0
    for r in range(0,runs):
        agent.initialiseOpinion(DecisionModel.LOGODDS, agents)
        if (agent.opinion > 0):
            success += 1
    
    print(success/runs)
    