'''
Created on 27 Jan 2017

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

import numpy as np
import scipy.special
import scipy.integrate as integrate
import math
import copy
import numpy.random as rand
import networkx as nx #@UnresolvedImport
import matplotlib.pyplot as plt
from DecNet.NetAgent import NetAgent
from DecNet.MyTypes import AgentType, DriftDistribution, NetworkType, DecisionModel, UpdateModel

class SpaceNet:
    def __init__(self):
        self.pos_dict = None
        self.orientations = None

class DecNet:
    
    DEBUG = False
    
    def __init__(self, netType, agentType, decModel, updateConf, numNodes, seed, debug=False):
        self.graph = []
        self.converged = False
        self.netType = netType
        self.agentType = agentType
        self.decModel = decModel
        self.updateConf = updateConf
        self.numNodes = numNodes
        self.agents = []
        self.randomSeed = seed*7
        self.DEBUG = debug
    
    def setNetworkParams(self, dynamic, args):
        if (self.netType == NetworkType.ERSOS_RENYI):
            self.linkProbability = args[0]
        elif(self.netType == NetworkType.BARABASI_ALBERT):
            self.numEdges = args[0]
        elif(self.netType == NetworkType.SPACE or self.netType == NetworkType.SOFT_RGG):
            self.areaSize = args[0]
            self.commRadius = args[1]
            self.periodicBounds = args[2]
        elif(self.netType == NetworkType.FROM_FILE or self.netType == NetworkType.FROM_FILE_FIXED_COMM):
            self.areaSize = args[0]
            self.coordFile = args[1]
            if(self.netType == NetworkType.FROM_FILE_FIXED_COMM):
                self.commRadius = args[2]
        self.dynamicNetwork = dynamic
        if self.dynamicNetwork:
            self.agentSpeed = args[-2]
            self.changeDirectionDegStddev = args[-1]
        
    def setSimpleAgent(self, accuracyMean, accuracyStdDev, truncatePoor):
        self.accuracyMean = accuracyMean
        self.accuracyStdDev = accuracyStdDev
        self.truncatePoor = truncatePoor
        
    def setDDMAgent(self, driftDistribution, baseDrift, noiseStdDev, interrogationTimeFromAccuracy, interrogationTime, prior, args):
        self.driftDistribution = driftDistribution
        self.baseDrift = baseDrift
        self.noiseStdDev = noiseStdDev
        self.DDMstart = 0
        self.dt = 0.01
        self.prior = prior
        if (self.driftDistribution == DriftDistribution.UNIFORM):
            self.randomDriftRangeMin = args[0]
            self.randomDriftRangeMax = args[1]
        elif (self.driftDistribution == DriftDistribution.NORMAL):
            self.driftStdDev = args[0]
        elif (self.driftDistribution == DriftDistribution.FROM_ACCURACY):
            self.accuracyMean = args[0]
            self.accuracyStdDev = args[1]
        # computing the interrogation time if that must give a given average accuracy as from Eq. (5.17) of Bogacz et al. Psy.Rev. 2006
        if interrogationTimeFromAccuracy:
            costMatrix = args[-1]
            expectedAcc = interrogationTime # I assume that the function parameter interrogationTime accepted the expected accuracy  
            self.interrogationTime = costMatrix[1]/costMatrix[0] * (2*expectedAcc - 1) * math.log( expectedAcc / (1-expectedAcc) ) / ( 2*math.log(expectedAcc / (1-expectedAcc)) - 1/expectedAcc + 1/(1-expectedAcc) )
            if (self.DEBUG): print("Computed interrogation time (for expected accuracy of " + str(expectedAcc) + ") is " + str(self.interrogationTime) )
        else:
            self.interrogationTime = interrogationTime
        
    
    def initNetwork(self):
        if (self.netType == NetworkType.FULLY_CONNECTED):
            self.graph = nx.complete_graph(self.numNodes) #np.repeat(0, self.numNodes)
        elif (self.netType == NetworkType.ERSOS_RENYI):
            self.graph = nx.erdos_renyi_graph(self.numNodes, self.linkProbability, self.randomSeed)
            i = 0
            while ( not nx.is_connected( self.graph ) ):
                if (self.DEBUG): print("Graph was not connected; Resampling!")
                i = i+1
                self.graph = nx.erdos_renyi_graph(self.numNodes, self.linkProbability, self.randomSeed*i*2211)
        elif (self.netType == NetworkType.BARABASI_ALBERT):
            self.graph = nx.barabasi_albert_graph(self.numNodes, self.numEdges, self.randomSeed)
#             self.graph = DriftDiffMod.Networks.powerlaw_cluster(self.numNodes, self.numEdges, 0, self.linkProbability, self.randomSeed)
        elif (self.netType == NetworkType.SPACE or self.netType == NetworkType.SOFT_RGG):
            connected = False
            while (not connected):
                self.spaceGraph = SpaceNet()
                self.spaceGraph.pos_dict = {}
                self.graph = nx.Graph()
                #g=nx.random_geometric_graph(self.numNodes, self.commRadius, pos={i: np.random.uniform(size=2) for i in range(self.numNodes) })
                for i in np.arange(0,self.numNodes):
                    self.spaceGraph.pos_dict[i] = np.array( [rand.uniform(0, self.areaSize), rand.uniform(0, self.areaSize)] )  
                    self.graph.add_node(i)
                    i += 1
                if self.dynamicNetwork:
                    connected = True
                    self.graph = nx.random_geometric_graph(self.numNodes, self.commRadius, pos=self.spaceGraph.pos_dict)
                    break # if the network is dynamic, the resampling for non-connected networks is not necessary
                else:
                    for agent in np.arange(0,self.numNodes):
                        for neigh, neighPos in self.spaceGraph.pos_dict.items():
                            if (self.periodicBounds):
                                if (not neigh == agent) and (self._distanceOnTorus(self.spaceGraph.pos_dict[agent][0], self.spaceGraph.pos_dict[agent][1], neighPos[0], neighPos[1]) <= self.commRadius):
                                    self.graph.add_edge(agent, neigh)
                            else:
                                if (not neigh == agent) and (np.linalg.norm(self.spaceGraph.pos_dict[agent]-neighPos) <= self.commRadius):
                                    self.graph.add_edge(agent, neigh)
                    
                    connected = nx.is_connected( self.graph )
                    if not connected:
                        if (self.DEBUG): print("Graph was not connected; Resampling!")
                    # else: self.spaceGraph.pos_dict=dict((i,self.spaceGraph.positions[i]) for i in range(0, len(self.spaceGraph.positions)))
                
        elif (self.netType == NetworkType.FROM_FILE_FIXED_COMM):
            file = open(self.coordFile, 'r')
            self.spaceGraph = SpaceNet()
            self.spaceGraph.pos_dict = {}
            line = file.readline()
            file.close()
            values = line.split('\t')[4:]
            self.numNodes = int(len(values)/4)
            #print(len(values))
            for idx in range(self.numNodes):
                self.spaceGraph.pos_dict[int(values[idx*4])] = [float(values[idx*4 +1])/self.areaSize, float(values[idx*4 +2])/self.areaSize]
            #print(self.spaceGraph.pos_dict)
            self.graph = nx.random_geometric_graph(self.numNodes, self.commRadius, pos=self.spaceGraph.pos_dict)
            if not nx.is_connected( self.graph ):
                if (self.DEBUG or not self.dynamicNetwork): print("WARNING! Input graph is not connected!")
        
        if (self.dynamicNetwork):
            self.spaceGraph.orientations = {}
            for n in self.getAllNodes():
                self.spaceGraph.orientations[n] = np.random.uniform(0,360)
    
    def getAllNodes(self):
        return self.graph.nodes()
#         if (self.netType == NetworkType.FULLY_CONNECTED):
#             return self.graph.nodes() #np.arange(self.numNodes)
#         elif (self.netType == NetworkType.ERSOS_RENYI):
#             return self.graph.nodes()
#         elif (self.netType == NetworkType.BARABASI_ALBERT):
#             return self.graph.nodes()
#         elif (self.netType == NetworkType.SPACE):
#             return self.graph.nodes()
        
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
    
    def initAgents(self, agentParams=[]):
        self.agents = []
        self.maxAcc = -float('inf')
        maxNeigh = -float('inf')
#         if (self.netType == NetworkType.FROM_FILE or self.netType == NetworkType.FROM_FILE_FIXED_COMM):
#             file = open(self.coordFile, 'r')
#             line = file.readline()
#             file.close()
#             values = line.split('\t')[4:]
#             for n in self.getAllNodes():
#                 args = []
#                 acc = float(values[n*4 +3])
#                 args.append( acc ) #accuracy
#                 self.agents.append( NetAgent(self.agentType, args, self.DEBUG) )
#                 if (self.agents[n].accuracy > self.maxAcc ):
#                     self.bestAccNode = n
#                     self.maxAcc  = self.agents[n].accuracy
#                 if (self.updateConf == UpdateModel.BELIEF_UP):
#                     self.agents[n].setBeliefEpsilon(beliefEpsilon)
#                     maxNeigh = max(maxNeigh, len(self.getNeighbours(n))) 
#         else:
        if (self.agentType == AgentType.SIMPLE):
            for n in self.getAllNodes():
                args = []
                acc = rand.normal(self.accuracyMean, self.accuracyStdDev)
                while (acc <= 0 or acc >=  1):
                    acc = rand.normal(self.accuracyMean, self.accuracyStdDev)
                if (self.truncatePoor and acc < 0.5):
                    acc = 1 - acc
                args.append( acc ) #accuracy
                self.agents.append( NetAgent(self.agentType, args, self.DEBUG) )
                if (self.agents[n].accuracy > self.maxAcc ):
                    self.bestAccNode = n
                    self.maxAcc  = self.agents[n].accuracy
                if (self.updateConf == UpdateModel.BELIEF_UP or self.updateConf == UpdateModel.FINITE_TIME):
                    if len(agentParams)>0:
                        self.agents[n].setBeliefEpsilon(agentParams[0])
                    if len(agentParams)>1:
                        self.agents[n].setFiniteTimeExponent(agentParams[1])
                    maxNeigh = max(maxNeigh, len(self.getNeighbours(n))) 
            if (self.DEBUG): print("MaxNeigh: " + str(maxNeigh) + " optimal Belief-Epsilon: " + str(1/maxNeigh if maxNeigh > 0 else np.Inf))
            
        elif (self.agentType == AgentType.DDM):
            # computing the optimal starting point, if prior probability is != 0.5
            if (self.driftDistribution == DriftDistribution.UNIFORM):
                self.driftStdDev = (self.randomDriftRangeMax - self.randomDriftRangeMin) / np.sqrt(12)
                meanDrift = self.baseDrift + ((self.randomDriftRangeMax + self.randomDriftRangeMin)/2.0)
            else:
                meanDrift = self.baseDrift
            if (not self.prior == 0.5):
                self.DDMstart = ( (self.noiseStdDev**2) / (2.0*meanDrift) ) * math.log(self.prior / (1-self.prior))
                print("Mean DDM is " + str(meanDrift) + ' and start value is ' + str(self.DDMstart))
                self.DDMstart = integrate.quad( self.integrationFuncStartDdmNormal, -np.inf, np.inf, args=( self.noiseStdDev, self.prior ) )[0]
                print('Integral of start values is ' + str(self.DDMstart))
                
            args = []
            args.append( 0 ) # space reserved for driftRate 
            args.append( self.noiseStdDev )
            args.append( self.DDMstart )
            args.append( self.dt )
            args.append( self.interrogationTime )
            args.append( self.prior )
                
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
                # revert negative drifts
                args[0] = abs(driftRate)
                self.agents.append( NetAgent(self.agentType, args, self.DEBUG) )
                
                if (self.driftDistribution  == DriftDistribution.FROM_ACCURACY):
                    self.agents[n].setPopMeanAndStdDev(self.accuracyMean, self.accuracyStdDev)
                else:
                    self.agents[n].setPopMeanAndStdDev(self.baseDrift, self.driftStdDev)
                
                # store the maximum drift
                if (self.agents[n].drift > self.maxAcc ):
                    self.bestAccNode = n
                    self.maxAcc  = self.agents[n].drift
          
    def initDecisions(self):
        if (self.netType == NetworkType.FROM_FILE or self.netType == NetworkType.FROM_FILE_FIXED_COMM):
            file = open(self.coordFile, 'r')
            line = file.readline()
            file.close()
            values = line.split('\t')[4:]
#             for n in self.getAllNodes():
#                 self.agents[n].opinion = int(values[n*4 +3])*2-1
                
        for n in self.getAllNodes(): 
            self.agents[n].initialiseOpinion(self.decModel, self.agents)
            if (self.netType == NetworkType.FROM_FILE or self.netType == NetworkType.FROM_FILE_FIXED_COMM): self.agents[n].opinion = int(values[n*4 +3])*2-1
        if (self.DEBUG and self.updateConf == UpdateModel.BELIEF_UP):
            expectedConsensus = 0
            for n in self.getAllNodes():
                expectedConsensus += (self.agents[n].opinion * self.agents[n].confidence)
            expectedConsensus /= self.numNodes
            print("Expected convergence to: " + str(expectedConsensus) )
    
    def getNeighbours(self, agent):
        if (self.netType == NetworkType.FULLY_CONNECTED):
            neighs = np.delete( np.arange(self.numNodes), agent)
        elif (self.netType == NetworkType.ERSOS_RENYI):
            neighs = list(nx.all_neighbors(self.graph, agent))
        elif (self.netType == NetworkType.BARABASI_ALBERT):
            neighs = list(nx.all_neighbors(self.graph, agent))
        elif (self.netType == NetworkType.SPACE):
            if (self.periodicBounds):
                neighs = []
                for neigh,neighPos in self.spaceGraph.pos_dict.items():
                    if (not neigh == agent) and (self._distanceOnTorus(self.spaceGraph.pos_dict[agent][0], self.spaceGraph.pos_dict[agent][1], neighPos[0], neighPos[1]) <= self.commRadius):
                        neighs.append(neigh)
                        self.graph.add_edge(agent, neigh)
            else:
                neighs = list(nx.all_neighbors(self.graph, agent))
        elif (self.netType == NetworkType.SOFT_RGG):
            neighs = []
            for neigh,neighPos in self.spaceGraph.pos_dict.items():
                if not neigh == agent:
                    if (self.periodicBounds):
                        dist = self._distanceOnTorus(self.spaceGraph.pos_dict[agent][0], self.spaceGraph.pos_dict[agent][1], neighPos[0], neighPos[1])
                    else:
                        dist = np.sqrt( (self.spaceGraph.pos_dict[agent][0] - neighPos[0])**2 + (self.spaceGraph.pos_dict[agent][1] - neighPos[1])**2) 
                        #np.linalg.norm(np.array(self.spaceGraph.pos_dict[agent][0], self.spaceGraph.pos_dict[agent][1])-np.array([neighPos[0], neighPos[1]]))
                    if (dist <= self.commRadius) and np.random.uniform() < (self.commRadius - dist)/self.commRadius:
                        neighs.append(neigh)
                        self.graph.add_edge(agent, neigh)
            
        elif (self.netType == NetworkType.FROM_FILE or self.netType == NetworkType.FROM_FILE_FIXED_COMM):
            neighs = list(nx.all_neighbors(self.graph, agent))
        return neighs        
    
    def _distanceOnTorus(self, x_1, y_1, x_2, y_2):
        """Returns the minimum distance calucalted on the torus given by periodic boundary conditions."""
        return np.sqrt(min(abs(x_1 - x_2), self.areaSize - abs(x_1 - x_2))**2 + 
                    min(abs(y_1 - y_2), self.areaSize - abs(y_1 - y_2))**2)
       
    def collectiveDecision(self, maxLoops, plot=None, plotTrajectory=False):
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
        
        if (self.DEBUG): logFile = open(self.logFilename, 'w')
        
        if plot is not None and plotTrajectory:
            self.trajectories = []
            for a in self.getAllNodes(): self.trajectories.append([ [self.spaceGraph.pos_dict[a][0],self.spaceGraph.pos_dict[a][1]] ])
             
        while (self.numNodes > self.countOpinion(1) and self.numNodes > self.countOpinion(-1) and count < maxLoops):
        #while (count < maxLoops):
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
                if (self.netType == NetworkType.SPACE or self.netType == NetworkType.SOFT_RGG or self.netType == NetworkType.FROM_FILE or self.netType == NetworkType.FROM_FILE_FIXED_COMM):
                    pos_layout = self.spaceGraph.pos_dict
                else:
                    pos_layout = nx.circular_layout(self.graph)
                nx.draw(self.graph, pos_layout, node_color=decisionColors, with_labels=True, ax=axs[0])
        
            # create a deep-copy of the current network (for synchronous communication)
            tmp_agents = copy.deepcopy(self.agents)
            if self.netType == NetworkType.SOFT_RGG: self.graph = nx.random_geometric_graph(self.numNodes, 0, pos=self.spaceGraph.pos_dict)
            # each node updates its opinion as the sum of its weighted opinion plus the neighbours' weighted opinions
            for n in self.getAllNodes():
                if (self.DEBUG): print('node ' + str(n))
                neighs = self.getNeighbours(n)
                if len(neighs) > 0:
                    self.agents[n].integrateInfoFromNeighbours(self.decModel, tmp_agents, neighs, self.updateConf )
                if (self.netType == NetworkType.FULLY_CONNECTED and self.updateConf != UpdateModel.BELIEF_UP): # if fully connected I stop after the first agent (the only WARNING is that I underestimate the conv time for majority-rand) 
                    #if (self.DEBUG): print('Collective choice is for option ' + str(self.agents[n].opinion))
#                     if (self.updateConf == UpdateModel.NO_UPDATE):
#                         self.agents[n].updateConfidenceOptimFull(tmp_agents, self.getNeighbours(n))
#                     estimatedAcc = self.agents[n].accuracy
                    if plot is not None:
                        plot.savefig()
                        plt.close(fig)
                    if (self.agents[n].opinion > 0):
                        return 1, len(self.agents), 0
                    else:
                        return 1, 0, len(self.agents)
        
            if (self.DEBUG): print("Positive nodes: " + str(self.countOpinion(1)) + " negative nodes: " + str(self.countOpinion(-1)) + " iteration: " + str(count) )
            
            if (self.DEBUG):
                line = str(count) + "\t" + str(self.countOpinion(1)/self.numNodes) + "\t" 
                for n in self.getAllNodes():
                    line += str(self.agents[n].confidence) + "\t"
                line += "\n"
                logFile.write(line)
            
            if plot is not None:
                if not plotTrajectory:
                    ## Draw final network
                    # update colors (red is opt+1, blue is opt-1) and plot
                    decisionColors = []
                    for n in self.getAllNodes():
                        if (self.agents[n].opinion > 0): decisionColors.append('r') 
                        elif  (self.agents[n].opinion < 0): decisionColors.append('b')
                        else: decisionColors.append('w')
                    nx.draw(self.graph, pos_layout, node_color=decisionColors, with_labels=True, ax=axs[1])
                else:
                    self.drawTrajectoryPlot(axs[1])
            
            # move agents and re-wire network 
            if self.dynamicNetwork: 
                for n in self.getAllNodes():
                    # add noise to orientation
                    self.spaceGraph.orientations[n] = (self.spaceGraph.orientations[n] + np.random.normal(0, self.changeDirectionDegStddev)) % 360
                    # update position 
                    currentPos = [self.spaceGraph.pos_dict[n][0],self.spaceGraph.pos_dict[n][1]]
                    currentPos[0] += self.agentSpeed * np.cos(np.deg2rad(self.spaceGraph.orientations[n])) 
                    currentPos[1] += self.agentSpeed * np.sin(np.deg2rad(self.spaceGraph.orientations[n]))
                    # out-of-bounds check
                    currentPos = self.solveOutOfBoundConditions(n, currentPos, plotTrajectory)
                    self.spaceGraph.pos_dict[n] = currentPos
                    if plotTrajectory: 
                        #if len(self.trajectories[n]) > 1: self.trajectories[n] = [ self.trajectories[n][-1] ]
                        #print("Appending to traj " + str(n) + " the current Pos: " + str(currentPos))
                        #print(self.trajectories)
                        self.trajectories[n].append( [currentPos[0],currentPos[1]] )
                # re-wire communication network
                self.graph = nx.random_geometric_graph(self.numNodes, self.commRadius, pos=self.spaceGraph.pos_dict)
            
            ## Print out the report
            if (self.DEBUG):
                for n in self.getAllNodes():
                    print('Node ' + str(n) + ' : ' + str(self.agents[n].opinion) + ' : ' + str(self.agents[n].confidence))
            
            if plot is not None:
                plot.savefig()
                plt.close(fig)
        
        if (self.DEBUG and self.updateConf == UpdateModel.BELIEF_UP):
            finalConfConsensus = 0
            for n in self.getAllNodes(): 
                finalConfConsensus += self.agents[n].confidence
                print("conf: " + str(self.agents[n].confidence) )
            finalConfConsensus /= self.numNodes
            print("Final convergence to: " + str(finalConfConsensus) )
        
        if (self.DEBUG): logFile.close()
        
        avg_degree = np.mean([ v for v   in dict(nx.degree(self.graph)).values() ]) if not self.netType == NetworkType.FULLY_CONNECTED else self.numNodes-1
        clust = nx.average_clustering(self.graph) if not self.netType == NetworkType.FULLY_CONNECTED else 1
        deg_stddev = np.std([ v for v   in dict(nx.degree(self.graph)).values() ]) if not self.netType == NetworkType.FULLY_CONNECTED else self.numNodes-1
        clust_stddev = np.std([ v for v   in dict(nx.clustering(self.graph)).values() ]) if not self.netType == NetworkType.FULLY_CONNECTED else self.numNodes-1
        confs = []
        for n in self.getAllNodes():
            confs.append(self.agents[n].confidence)
        return count, self.countOpinion(1), self.countOpinion(-1), avg_degree, clust, deg_stddev, clust_stddev, np.mean(confs), np.std(confs)
        
    def solveOutOfBoundConditions(self, agent, currentPos, saveTrajectory):
        if self.periodicBounds:
            if currentPos[0] > self.areaSize: currentPos[0] -= self.areaSize
            if currentPos[0] < 0: currentPos[0] += self.areaSize
            if currentPos[1] > self.areaSize: currentPos[1] -= self.areaSize
            if currentPos[1] < 0: currentPos[1] += self.areaSize
        else:
            if currentPos[0] >= self.areaSize or currentPos[0] <= 0 or currentPos[1] >= self.areaSize or currentPos[1] <= 0:
                collision = True
                travel_left = self.agentSpeed
                prevPos = _point2D(self.spaceGraph.pos_dict[agent][0],self.spaceGraph.pos_dict[agent][1])
                southWestPoint = _point2D(0,0)
                southEastPoint = _point2D(self.areaSize,0)
                northWestPoint = _point2D(0,self.areaSize)
                northEastPoint = _point2D(self.areaSize,self.areaSize)
                loops = 0
                while collision:
                    loops += 1
                    nextPos = _point2D(currentPos[0],currentPos[1])
                    southColl = _doIntersect(southWestPoint,southEastPoint,prevPos,nextPos) and prevPos.y > 0 # southWall
                    eastColl = _doIntersect(southEastPoint,northEastPoint,prevPos,nextPos) and prevPos.x < self.areaSize # eastWall
                    northColl = _doIntersect(northWestPoint,northEastPoint,prevPos,nextPos) and prevPos.y < self.areaSize # northWall
                    westColl = _doIntersect(southWestPoint,northWestPoint,prevPos,nextPos) and prevPos.x > 0 # westWall
                    if southColl or eastColl or northColl or westColl:
                        tmp_slope = (nextPos.y-prevPos.y)/(nextPos.x-prevPos.x)
                        tmp_coeff = nextPos.y - (tmp_slope * nextPos.x)
                        if southColl:
                            currentPos[0] = -tmp_coeff/tmp_slope
                            currentPos[1] = 0
                            self.spaceGraph.orientations[agent] = 360 - self.spaceGraph.orientations[agent]
                        elif northColl:
                            currentPos[0] = (self.areaSize-tmp_coeff)/tmp_slope
                            currentPos[1] = self.areaSize
                            self.spaceGraph.orientations[agent] = 360 - self.spaceGraph.orientations[agent]
                        elif westColl:
                            currentPos[0] = 0
                            currentPos[1] = tmp_coeff
                            self.spaceGraph.orientations[agent] = 180 - self.spaceGraph.orientations[agent]
                        elif eastColl:
                            currentPos[0] = self.areaSize
                            currentPos[1] = tmp_slope*self.areaSize + tmp_coeff
                            self.spaceGraph.orientations[agent] = 180 - self.spaceGraph.orientations[agent]
                        travelled = np.linalg.norm(np.array([prevPos.x, prevPos.y])-np.array(currentPos))
                        if saveTrajectory: self.trajectories[agent].append( [currentPos[0],currentPos[1]] )
                        #print("* * * COLLISION * * * " + str(currentPos))
                        # bounce back
                        travel_left -= travelled
                        prevPos = _point2D(currentPos[0],currentPos[1])
                        currentPos[0] += travel_left * np.cos(np.deg2rad(self.spaceGraph.orientations[agent])) 
                        currentPos[1] += travel_left * np.sin(np.deg2rad(self.spaceGraph.orientations[agent]))
                        if (loops > 10): self.spaceGraph.orientations[agent] += np.random.normal(0, 1)
                    else:
                        collision = False
        return currentPos
        
    def drawTrajectoryPlot(self, ax):
        show_links = True
        show_radius = True
        show_traj = True
        
        if show_traj:
            #print(self.trajectories)
            for j in self.getAllNodes():
                xs = [x[0] for x in self.trajectories[j]]
                ys = [x[1] for x in self.trajectories[j]]
                #print(xs)
                ax.plot(xs,ys)
            
        xs = [x[0] for x in sorted(self.spaceGraph.pos_dict.values(), key=str)]
        ys = [x[1] for x in sorted(self.spaceGraph.pos_dict.values(), key=str)]
        ax.scatter(xs, ys)
        #ax = plt.gca()
#                    ax.cla() # clear things for fresh plot
        for a,pos in self.spaceGraph.pos_dict.items():
            if show_links:
                for n in self.getNeighbours(a):
                    ax.plot( [pos[0], self.spaceGraph.pos_dict[n][0]], [pos[1], self.spaceGraph.pos_dict[n][1]])
            if show_radius:
                if not self.getNeighbours(a):
                    comm_circle = plt.Circle((pos[0], pos[1]), self.commRadius, color='r', fill=False)
                else:
                    comm_circle = plt.Circle((pos[0], pos[1]), self.commRadius, color='b', fill=False)
                ax.add_artist(comm_circle)
        ax.axis([0, self.areaSize, 0, self.areaSize])
        ax.set_aspect('equal')
        #plt.show()
        
def _onSegment(p,q,r):
    if (q.x<= max(p.x,r.x) and q.x >= min(p.x,r.x) and q.y<= max(p.y,r.y) and q.y>=min(p.y,r.y)):
        return True
    return False

def _orientation(p,q,r):
    """Find the orientation of ordered triplet (p, q, r) 
    Returns: 0 if p, q, and r are collinear; 1 for Clockwise; and 2 for Counterclockwise"""
    val = (q.y - p.y)* (r.x - q.x) - (q.x-p.x)*(r.y-q.y)
    if val==0:
        return 0
    if val>0:
        return 1
    else:
        return 2
    
def _doIntersect(p1,q1,p2,q2):
    """check intersection between segment1 (p1,q1) and segment2 (p2, q2)
    Code adapted from: https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect"""
    o1 = _orientation(p1,q1,p2)
    o2 = _orientation(p1,q1,q2)
    o3 = _orientation(p2,q2,p1)
    o4 = _orientation(p2,q2,q1)
    
    if (not o1 == o2) and (not o3 == o4):
        return True
    if (o1 == 0 and _onSegment(p1,p2,q1)):
        return True
    if (o2 == 0 and _onSegment(p1,q2,q1)):
        return True
    if (o3 == 0 and _onSegment(p2,p1,q2)):
        return True
    if (o4 == 0 and _onSegment(p2,q1,q2)):
        return True
    return False

class _point2D:
    def __init__(self,v1,v2):
        self.x = v1
        self.y = v2



