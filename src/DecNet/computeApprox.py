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
    accMean = 0.6
    accStdev = 0.12
    noiseStdDev = 0.5
    interrogationTime = 2
    DDMstart = 0
    dt = 0.01
    seed = 221189
    runs = 10000
    
    rand.seed(seed)
    driftRate = computeDriftFromAccuracy(accMean, noiseStdDev, interrogationTime)
    
    args = []
    args.append( driftRate ) 
    args.append( noiseStdDev )
    args.append( DDMstart )
    args.append( dt )
    args.append( interrogationTime )
     
    agent = NetAgent(AgentType.DDM, args, True)
    
    logodds = []
    agent.setMeanAccuracyAndStdDev(accMean, accStdev)
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
    
    plt.ylabel('Confidence $c^{dL}$')
    plt.xlabel('Integrated value $x_i$')
    plt.plot(xrange, logodds, 'ro', ms=10.0)
    plt.plot(xrange, fittedline(xrange), 'k--', linewidth=3.0, label="Linear")
    plt.plot(xrange, fittedline2(xrange), 'b', linewidth=3.0, label="Quadratic")
    #plt.axis([0, 6, 0, 20])
#     blue_line = mlines.Line2D([], [], color='blue', marker='*',
#                           markersize=15, label='Blue stars')
#     plt.legend(handles=[blue_line])
    plt.legend(loc=(0.1, 0.75), borderaxespad=0.)
    plt.savefig('~/Google Drive/DiODe/Manuscripts/DDM-on-Net/plots/approxCurves.pdf')
    plt.show()
    
    exit()
    