'''
Created on 27 Jan 2017

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

import numpy as np
import scipy.special
import math
import copy
import numpy.random as rand
import networkx as nx #@UnresolvedImport
import matplotlib.pyplot as plt
from DecNet.NetAgent import NetAgent
from DecNet.MyTypes import NetworkType
from DecNet.MyTypes import DecisionModel
from DecNet.MyTypes import DriftDistribution
from DecNet.MyTypes import AgentType

class DecNet:
    
    DEBUG = False
    
    def __init__(self, netType, agentType, decModel, numNodes, seed, debug=False):
        self.graph = []
        self.converged = False
        self.netType = netType
        self.agentType = agentType
        self.decModel = decModel
        self.numNodes = numNodes
        self.agents = []
        self.randomSeed = seed*7
        self.DEBUG = debug
    
    def setNetworkParams(self, args):
        if (self.netType == NetworkType.ERSOS_RENYI):
            self.linkProbability = args[0]
        elif(self.netType == NetworkType.BARABASI_ALBERT):
            self.numEdges = args[0]
        
    def setSimpleAgent(self, accuracyMean, accuracyStdDev):
        self.accuracyMean = accuracyMean
        self.accuracyStdDev = accuracyStdDev
        
    def setDDMAgent(self, driftDistribution, baseDrift, noiseStdDev, interrogationTime, args):
        self.driftDistribution = driftDistribution
        self.baseDrift = baseDrift
        self.noiseStdDev = noiseStdDev
        self.DDMstart = 0
        self.dt = 0.01
        self.interrogationTime = interrogationTime
        if (self.driftDistribution == DriftDistribution.UNIFORM):
            self.randomDriftRangeMin = args[0]
            self.randomDriftRangeMax = args[1]
        elif (self.driftDistribution == DriftDistribution.NORMAL):
            self.driftStdDev = args[0]
        elif (self.driftDistribution == DriftDistribution.FROM_ACCURACY):
            self.accuracyMean = args[0]
            self.accuracyStdDev = args[1]
        
    
    def initNetwork(self):
        if (self.netType == NetworkType.FULLY_CONNECTED):
            self.graph = nx.complete_graph(self.numNodes) #np.repeat(0, self.numNodes)
        elif (self.netType == NetworkType.ERSOS_RENYI):
            self.graph = nx.erdos_renyi_graph(self.numNodes, self.linkProbability, self.randomSeed)
        elif (self.netType == NetworkType.BARABASI_ALBERT):
            self.graph = nx.barabasi_albert_graph(self.numNodes, self.numEdges)
#             self.graph = DriftDiffMod.Networks.powerlaw_cluster(self.numNodes, self.numEdges, 0, self.linkProbability, self.randomSeed)
    
    def getAllNodes(self):
        if (self.netType == NetworkType.FULLY_CONNECTED):
            return self.graph.nodes() #np.arange(self.numNodes)
        elif (self.netType == NetworkType.ERSOS_RENYI):
            return self.graph.nodes()
        elif (self.netType == NetworkType.BARABASI_ALBERT):
            return self.graph.nodes()
        
    def countDecided(self):
        count = 0
        for a in self.agents:
            if (a.opinion != 0):
                count += 1
        return count

    def countOpinion(self, op):
        count = 0
        for a in self.agents:
            if (a.opinion == op):
                count += 1
        return count
    
    def computeDriftFromAccuracy(self, accuracy, noise, interrogationTime):
        return math.sqrt(2) * noise * scipy.special.erfcinv(2 - 2*accuracy)/math.sqrt(interrogationTime)
    
    def initAgents(self):
        self.agents = []
        self.maxAcc = -float('inf')
        if (self.agentType == AgentType.SIMPLE):
            for n in self.getAllNodes():
                args = []
                acc = rand.normal(self.accuracyMean, self.accuracyStdDev)
                while (acc < 0 or acc >=  1):
                    acc = rand.normal(self.accuracyMean, self.accuracyStdDev)
                args.append( acc ) #accuracy
                self.agents.append( NetAgent(self.agentType, args, self.DEBUG) )
                if (self.agents[n].accuracy > self.maxAcc ):
                    self.bestAccNode = n
                    self.maxAcc  = self.agents[n].accuracy 
        elif (self.agentType == AgentType.DDM):
            args = []
            args.append( 0 ) # space reserved for driftRate 
            args.append( self.noiseStdDev )
            args.append( self.DDMstart )
            args.append( self.dt )
            args.append( self.interrogationTime )
                
            for n in self.getAllNodes():
                # randomly select a drift
                if (self.driftDistribution == DriftDistribution.UNIFORM):
                    driftRate = self.baseDrift + rand.uniform(self.randomDriftRangeMin, self.randomDriftRangeMax)
                elif (self.driftDistribution == DriftDistribution.NORMAL):
                    driftRate = rand.normal(self.baseDrift, self.driftStdDev)
                elif (self.driftDistribution == DriftDistribution.FROM_ACCURACY):
                    acc = rand.normal(self.accuracyMean, self.accuracyStdDev)
                    while (acc < 0 or acc >=  1):
                        acc = rand.normal(self.accuracyMean, self.accuracyStdDev)
                    driftRate = self.computeDriftFromAccuracy(acc, self.noiseStdDev, self.interrogationTime)
                
                args[0] = driftRate
                self.agents.append( NetAgent(self.agentType, args, self.DEBUG) )
                
                # store the maximum drift
                if (self.agents[n].drift > self.maxAcc ):
                    self.bestAccNode = n
                    self.maxAcc  = self.agents[n].drift
            
          
    def initDecisions(self):
        for n in self.getAllNodes():
            self.agents[n].initialiseOpinion(self.decModel, self.agents)
               
    
    def getNeighbours(self, agent):
        if (self.netType == NetworkType.FULLY_CONNECTED):
            neighs = np.delete( np.arange(self.numNodes), agent)
        elif (self.netType == NetworkType.ERSOS_RENYI):
            neighs = list(nx.all_neighbors(self.graph, agent))
        elif (self.netType == NetworkType.BARABASI_ALBERT):
            neighs = list(nx.all_neighbors(self.graph, agent))
        return neighs
    
            
    def collectiveDecision(self, maxLoops, plot=None):
        if (self.decModel == DecisionModel.BEST_ACC):
            if (self.agents[self.bestAccNode].opinion == 1):
                allPos = len(self.agents) 
            else:
                allPos = 0
            allNeg = len(self.agents) - allPos
            return 1, allPos, allNeg
        
        count = 0
        ## Print out the report
        if (self.DEBUG):
            for n in self.getAllNodes():
                print('Node ' + str(n) + ' : ' + str(self.agents[n].opinion) + ' : ' + str(self.agents[n].confidence))
                
        while (self.numNodes > self.countOpinion(1) and self.numNodes > self.countOpinion(-1) and count < maxLoops):
            count += 1
            if plot is not None: 
                ## Draw initial network
                fig, axs = plt.subplots(1,2)
                axs[0].set_title('Iteration n.' + str(count))
                decisionColors=[]
                for n in self.getAllNodes():
                    if (self.agents[n].opinion > 0): decisionColors.append('r') 
                    elif  (self.agents[n].opinion < 0): decisionColors.append('b')
                    else: decisionColors.append('w')
                nx.draw(self.graph, nx.circular_layout(self.graph), node_color=decisionColors, with_labels=True, ax=axs[0])
        
            # create a deep-copy of the current network (for synchronous communication)
            tmp_agents = copy.deepcopy(self.agents)
            # each node updates its opinion as the sum of its weighted opinion plus the neighbours' weighted opinions
            for n in self.getAllNodes():
                if (self.DEBUG): print('node ' + str(n))
                self.agents[n].integrateInfoFromNeighbours(self.decModel, tmp_agents, self.getNeighbours(n) )
                if (self.DEBUG): print('My new "confidence" (node:' + str(n) + ') is ' + str(self.agents[n].confidence))
                if (self.netType == NetworkType.FULLY_CONNECTED): # if fully connected I stop after the first agent (the only WARNING is that I underestimate the conv time for majority-rand) 
                    #if (self.DEBUG): print('Collective choice is for option ' + str(self.agents[n].opinion))
                    if (self.agents[n].opinion > 0):
                        return 1, len(self.agents), 0 
                    else:
                        return 1, 0, len(self.agents)
        
            if (self.DEBUG): print("Positive nodes: " + str(self.countOpinion(1)) + " negative nodes: " + str(self.countOpinion(-1)) + " iteration: " + str(count) )
        
            if plot is not None:
                ## Draw final network
                # update colors (red is opt+1, blue is opt-1) and plot
                decisionColors = []
                for n in self.getAllNodes():
                    if (self.agents[n].opinion > 0): decisionColors.append('r') 
                    elif  (self.agents[n].opinion < 0): decisionColors.append('b')
                    else: decisionColors.append('w')
                nx.draw(self.graph, nx.circular_layout(self.graph), node_color=decisionColors, with_labels=True, ax=axs[1])
        
            ## Print out the report
            if (self.DEBUG):
                for n in self.getAllNodes():
                    print('Node ' + str(n) + ' : ' + str(self.agents[n].opinion) + ' : ' + str(self.agents[n].confidence))
            
            if plot is not None:
                plot.savefig()
                plt.close(fig)
        return count, self.countOpinion(1), self.countOpinion(-1)
        
        