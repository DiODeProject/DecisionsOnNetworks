'''
Created on 09 Mar 2017

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

import scipy.special
import math
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
from DecNet.NetAgent import NetAgent
from DecNet.MyTypes import AgentType 
from DecNet.MyTypes import DecisionModel 

def computeDriftFromAccuracy(accuracy, noise, interrogationTime):
        return math.sqrt(2) * noise * scipy.special.erfcinv(2 - 2*accuracy)/math.sqrt(interrogationTime)

if __name__ == '__main__':
    accuracy = 0.8217
    noiseStdDev = 0.5
    interrogationTime = 2
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
    
    pdf = lambda decVar, drift, noise, time: np.exp( (-(decVar - drift * time)**2) / (2 * time * noise * noise) ) / math.sqrt(2 * math.pi * time * noise * noise)
    print( pdf(0.4, 0.1, noiseStdDev, interrogationTime) )
    
    logodds = []
    agent.setMeanAccuracyAndStdDev(0.6, 0.12)
    xrange = np.arange(0, 4, 0.1)
    for x in xrange:
        logodds.append( agent.logOddsDistribution(x) )
    print(logodds)
    
    fit = np.polyfit(xrange, logodds, 1)
    print(np.poly1d(fit))
    fittedline = np.poly1d(fit)
    fit2 = np.polyfit(xrange, logodds, 2)
    print(np.poly1d(fit2))
    fittedline2 = np.poly1d(fit2)
    #fittedline3 = np.poly1d([0.1663, 0.5309, 0.1238])
    
    plt.plot(xrange, logodds, 'ro')
    plt.plot(xrange, fittedline(xrange), 'k')
    plt.plot(xrange, fittedline2(xrange), 'b')
    #plt.axis([0, 6, 0, 20])
    plt.show()
    
    exit()
    
    success = 0
    for r in range(0,runs):
        agent.initialiseOpinion(DecisionModel.LOGODDS_PERFECT, agents)
        if (agent.opinion > 0):
            success += 1
    
    print(success/runs)
    