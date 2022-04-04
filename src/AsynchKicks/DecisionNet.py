'''
Created on 15 Mar 2019

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

import numpy as np
import scipy.special
import scipy.optimize
from scipy.stats import norm
import scipy.integrate as integrate
import math
import numpy.random as rand
import networkx as nx #@UnresolvedImport
import matplotlib.pyplot as plt
from AsynchKicks.DdmAgent import DdmAgent
from AsynchKicks.MyTypes import AgentType, DriftDistribution, NetworkType, UpdateModel

class SpaceNet:
    def __init__(self):
        self.pos_dict = None
        self.orientations = None
        
class CascadeLog:
    def __init__(self, rank, cascadeSize, opinion):
        self.rank = rank
        self.cascadeSize = cascadeSize
        self.opinion = opinion
        
    def __str__(self): 
        return str("rank:" + str(self.rank) + ", casc:" + str(self.cascadeSize) + ", op:" + str(self.opinion))
    
    def __repr__(self):
        return self.__str__() 

class DecNet:
    
    DEBUG = False
    
    def __init__(self, netType, agentType, updateModel, numNodes, seed, debug=False):
        self.graph = []
        self.converged = False
        self.netType = netType
        self.agentType = agentType
        self.updateModel = updateModel
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
        
    def setDDMAgent(self, driftDistribution, baseDrift, noiseStdDev, dt, threshold, prior, useBayesRisk, misinformed, args):
        self.driftDistribution = driftDistribution
        self.baseDrift = baseDrift
        self.noiseStdDev = noiseStdDev
        self.DDMstart = 0
        self.dt = dt
        self.prior = prior
        self.misinformed = misinformed
        if (self.driftDistribution == DriftDistribution.UNIFORM):
            self.randomDriftRangeMin = args[0]
            self.randomDriftRangeMax = args[1]
            meanDrift = self.baseDrift + ((self.randomDriftRangeMax + self.randomDriftRangeMin)/2.0)
        elif (self.driftDistribution == DriftDistribution.NORMAL):
            self.driftStdDev = args[0]
            meanDrift = self.baseDrift
        elif (self.driftDistribution == DriftDistribution.FROM_ACCURACY):
            self.accuracyMean = args[0]
            self.accuracyStdDev = args[1]
        if useBayesRisk:
            self.costMatrix = args[-1]
            integrateOverAllValues = True
            if integrateOverAllValues:
                integrationLowerLimit = -np.inf if self.misinformed else 0
                normalisationOfNormProbForResampling = 1 if self.misinformed else integrate.quad( norm.pdf, 0, np.inf, args=(self.baseDrift, self.driftStdDev))[0]
                integrResult = integrate.quad( self.integrationFuncDdmThreshNormal, integrationLowerLimit, np.inf, args=( self.costMatrix, self.noiseStdDev, normalisationOfNormProbForResampling ) )
                self.threshold = integrResult[0]
                if (self.DEBUG): print("Computed (with integral with error: " + str(integrResult[1]) + ") threshold with BR cost matrix is: " + str(self.threshold) )
            else:
                def compute_BR_opt_thresh(thresh): # from Eq. (5.6) of Bogacz et al. Psy.Rev. 2006
                    return self.costMatrix[1]/self.costMatrix[0] * 2 * (meanDrift**2) / (self.noiseStdDev**2) - 4 * meanDrift * thresh / (self.noiseStdDev**2) + np.exp( - 2 * meanDrift * thresh / (self.noiseStdDev**2) ) -  np.exp( 2 * meanDrift * thresh / (self.noiseStdDev**2) )            
                self.threshold = scipy.optimize.fsolve( compute_BR_opt_thresh, 0.25*meanDrift*self.costMatrix[1]/self.costMatrix[0], maxfev=1000 )[0] # second parameter is the starting point which is set to the approximation of Eq. (5.7) of Bogacz et al. Psy.Rev. 2006 
                if (self.DEBUG): print("Computed threshold with BR cost matrix is: " + str(self.threshold) )
                #if (self.DEBUG): print("CHECK: " + str( compute_BR_opt_thresh(self.threshold)) )
        else:
            self.threshold = threshold
        
    
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
    
    def integrationFuncDdmThreshNormal(self, drift, costMatrix, noiseStdDev, normalisationOfNormProbForResampling):
        prob = norm.pdf(drift, self.baseDrift, self.driftStdDev) / normalisationOfNormProbForResampling
        #if not misinformed: drift = abs(drift)
        def compute_BR_opt_thresh(thresh): # from Eq. (5.6) of Bogacz et al. Psy.Rev. 2016
            return costMatrix[1]/costMatrix[0] * 2 * (drift**2) / (noiseStdDev**2) - 4 * drift * thresh / (noiseStdDev**2) + np.exp( - 2 * drift * thresh / (noiseStdDev**2) ) -  np.exp( 2 * drift * thresh / (noiseStdDev**2) )             
        pdf = scipy.optimize.fsolve( compute_BR_opt_thresh, 0.25*drift*costMatrix[1]/costMatrix[0], maxfev=5000 )[0] # second parameter is the starting point which is set to the approximation of Eq. (5.7) of Bogacz et al. Psy.Rev. 2016
        return pdf*prob
    
    def integrationFuncStartDdmNormal(self, drift, noiseStdDev, prior, normalisationOfNormProbForResampling):
        prob = norm.pdf(drift, self.baseDrift, self.driftStdDev) / normalisationOfNormProbForResampling
#         if misinformed:
#             pdf = ( (noiseStdDev**2) / (2.0*drift) ) * math.log(prior / (1-prior)) # possible negative drift 
#         else:  
#             pdf = ( (noiseStdDev**2) / (2.0*abs(drift)) ) * math.log(prior / (1-prior)) # assuming always positive drifts (flipped if negative)
        pdf = ( (noiseStdDev**2) / (2.0*drift) ) * math.log(prior / (1-prior))
        return pdf*prob
    
    def initAgents(self):
        self.agents = []
        self.maxAcc = -float('inf')
        if (self.agentType == AgentType.SIMPLE):
            for n in self.getAllNodes():
                args = []
                acc = rand.normal(self.accuracyMean, self.accuracyStdDev)
                while (acc <= 0 or acc >=  1):
                    acc = rand.normal(self.accuracyMean, self.accuracyStdDev)
                if (self.truncatePoor and acc < 0.5):
                    acc = 1 - acc
                args.append( acc ) #accuracy
                self.agents.append( DdmAgent(n, self.agentType, self.updateModel, args, self.DEBUG) )
                if (self.agents[n].accuracy > self.maxAcc ):
                    self.bestAccNode = n
                    self.maxAcc  = self.agents[n].accuracy
            
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
                integrationLowerLimit = -np.inf if self.misinformed else 0
                normalisationOfNormProbForResampling = 1 if self.misinformed else integrate.quad( norm.pdf, 0, np.inf, args=(self.baseDrift, self.driftStdDev))[0]
                self.DDMstart = integrate.quad( self.integrationFuncStartDdmNormal, integrationLowerLimit, np.inf, args=( self.noiseStdDev, self.prior, normalisationOfNormProbForResampling ) )[0]
                print('Integral of start values is ' + str(self.DDMstart))
                
            args = []
            args.append( 0 ) # space reserved for driftRate 
            args.append( self.noiseStdDev )
            args.append( self.DDMstart )
            args.append( self.dt )
            args.append( self.threshold )
            args.append( self.prior )
            args.append( self.misinformed )
                
            for n in self.getAllNodes():
                # randomly select a drift
                if (self.driftDistribution == DriftDistribution.UNIFORM):
                    driftRate = self.baseDrift + rand.uniform(self.randomDriftRangeMin, self.randomDriftRangeMax)
                elif (self.driftDistribution == DriftDistribution.NORMAL):
                    driftRate = rand.normal(self.baseDrift, self.driftStdDev)
                    # if misinformed, resample the drift until it is positive
                    if not self.misinformed:
                        while driftRate < 0:
                            driftRate = rand.normal(self.baseDrift, self.driftStdDev)
                elif (self.driftDistribution == DriftDistribution.FROM_ACCURACY):
                    acc = rand.normal(self.accuracyMean, self.accuracyStdDev)
                    while (acc < 0 or acc >=  1):
                        acc = rand.normal(self.accuracyMean, self.accuracyStdDev)
                    driftRate = self.computeDriftFromAccuracy(acc, self.noiseStdDev, self.interrogationTime)
                ## revert negative drifts
                #args[0] = driftRate if self.misinformed else abs(driftRate)
                args[0] = driftRate  
                self.agents.append( DdmAgent(n, self.agentType, self.updateModel, args, self.DEBUG) )
                
                if (self.driftDistribution  == DriftDistribution.FROM_ACCURACY):
                    self.agents[n].setMeanAndStdDev(self.accuracyMean, self.accuracyStdDev)
                else:
                    self.agents[n].setMeanAndStdDev(self.baseDrift, self.driftStdDev)
                
                # store the maximum drift
                if (self.agents[n].drift > self.maxAcc ):
                    self.bestAccNode = n
                    self.maxAcc  = self.agents[n].drift
            tmp_sorting = {}
            for a,agent in enumerate(self.agents):
                tmp_sorting[a] = agent.drift
            self.sortedAgentsByDescDrift = sorted(tmp_sorting.items(), key=lambda val: abs(val[1]))
            self.sortedAgentsByDescDrift.reverse()
            if (self.DEBUG): print("Sorted list of drifts: " +str(self.sortedAgentsByDescDrift))
    
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
    
    def plotNet(self, plot, time):
        fig, axs = plt.subplots(1,1)
        axs.set_title('Time ' + str(time))
        decisionColors=[]
        for n in self.getAllNodes():
            if (self.agents[n].opinion > 0): decisionColors.append('r') 
            elif  (self.agents[n].opinion < 0): decisionColors.append('g')
            else: decisionColors.append('y')
        if (self.netType == NetworkType.SPACE or self.netType == NetworkType.SOFT_RGG or self.netType == NetworkType.FROM_FILE or self.netType == NetworkType.FROM_FILE_FIXED_COMM):
            pos_layout = self.spaceGraph.pos_dict
        else:
            pos_layout = nx.circular_layout(self.graph)
        nx.draw(self.graph, pos_layout, node_color=decisionColors, with_labels=True, ax=axs)
        plot.savefig()
        plt.close(fig)
        
    def plotDDM(self, plot, _time, aID):
        fig, axs = plt.subplots(1,1)
        axs.set_title('Agent ' + str(aID) + "[drift: " + str(self.agents[aID].drift) +"]")
        axs.axhline(y=self.threshold,  xmin=0, xmax=1, ls='dashed', color='r')
        axs.axhline(y=0, xmin=0, xmax=1, ls='dashed', color='lightgray')
        axs.axhline(y=-self.threshold, xmin=0, xmax=1, ls='dashed', color='g')
        axs.set_ylim(-self.threshold*1.05, self.threshold*1.05)
        plt.plot( self.agents[aID].DDMintegration_time, self.agents[aID].DDMintegration, color='b')
        plot.savefig()
        plt.close(fig)
        
    def logNewKickSize(self, kickerOpinion, time):
        if (self.updateModel == UpdateModel.NO_UPDATE):
            kicksize = 0
        if (self.updateModel == UpdateModel.THRESH_KICK):
            kicksize = self.threshold
        if (self.updateModel == UpdateModel.CONF_KICK):
            priorComponent = 0 if (self.prior == 0.5) else math.log( self.prior/ (1-self.prior))  # if the sing of kickerOpinion is negative, through math.copysign I change the sing of the priorComponent (which is equivalent to power to -1 the log argument)
            kicksize = kickerOpinion * self.agents[0].estimateConfFromDistribution(self.threshold, time) - math.copysign(priorComponent, kickerOpinion)
        return kicksize
    
    def collectiveDecision(self, maxTime, plot=None):
        time = 0
        for agent in self.agents: agent.initialiseDDM() # initialise the DDM
        if plot is not None: self.plotNet(plot, time) ## Draw initial network
        
        countDecided = self.countDecided() # using this counter rather than function countDecided() to save time 
        kickSizes = []
        cascadeCounts = []
        cascadeLog = []
        sumOfIndividualBayesRisks = 0
        while (time < maxTime):
            time += self.dt
            agentsWhoKicked = []
            # integrate individual evidence via DDM
            for a,agent in enumerate(self.agents):
                if agent.opinion == 0: # decided agents do not integrate evidence
                    if not agent.integrateIndividualEvidence(time) == 0: # if not 0, the agent kicked
                        countDecided += 1
                        agentsWhoKicked.append(a) 
                        sumOfIndividualBayesRisks += (agent.opinion==-1)*self.costMatrix[1] + time*self.costMatrix[0] 
                        kickSizes.append(self.logNewKickSize(agent.opinion, time))
                        if (self.DEBUG): print('Node ' + str(a) + ' kicked at time ' + str(time) + ' to ' + str(agent.opinion))
                        if plot is not None: self.plotDDM(plot, time, a)
            if plot is not None and len(agentsWhoKicked) > 0: self.plotNet(plot, time)
            if len(agentsWhoKicked) > 0:
                for independentKicker in agentsWhoKicked:
                    ikPosition = 0
                    for a,sortedAgent in enumerate(self.sortedAgentsByDescDrift):
                        if independentKicker == sortedAgent[0]:
                            ikPosition += a+1
                    ikOpinion = self.agents[independentKicker].opinion
                ikPosition = ikPosition/len(agentsWhoKicked)
                if ikPosition <= 0: print (" * * * * * * ERROR IN RETRIEVING THE AGENT POSTION * * * * * *"); return
                elif (self.DEBUG): print("Kicker(s) position is: " + str(ikPosition))
            cascadeCount = len(agentsWhoKicked)
            while len(agentsWhoKicked) > 0:
                np.random.shuffle(agentsWhoKicked) # todo: handle the case of two kicks in opposite directions at the same time... shall they be processed sequentially or synchronously?
                cascadeKicks = []
                for kicker in agentsWhoKicked: # information is sent to neighbours
                    neighs = self.getNeighbours(kicker)
                    for neigh in neighs:
                        if self.agents[neigh].opinion == 0: # decided agents do not receive kicks
                            if not self.agents[neigh].kick(self.agents[kicker].opinion, time) == 0: # if not 0, the agent kicked
                                countDecided += 1
                                cascadeKicks.append(neigh) # a kick can create a cascade of kicks
                                sumOfIndividualBayesRisks += (self.agents[neigh]==-1)*self.costMatrix[1] + time*self.costMatrix[0]
                                if (self.DEBUG): print('CASCADE: Node ' + str(neigh) + ' kicked at time ' + str(time) + ' to ' + str(self.agents[neigh].opinion))
                            if plot is not None: self.plotDDM(plot, time, neigh)
                agentsWhoKicked.clear()
                agentsWhoKicked.extend(cascadeKicks)
                cascadeCount += len(agentsWhoKicked)
                if plot is not None: self.plotNet(plot, time)
            #if (self.numNodes == self.countOpinion(1) or self.numNodes == self.countOpinion(-1)): # consensus reached
            if cascadeCount > 0: 
                cascadeCounts.append(cascadeCount)
                if (self.DEBUG): print("Cascade: " + str(cascadeCounts))
                cascadeLog.append(CascadeLog(ikPosition, cascadeCount, ikOpinion))
            if self.numNodes == countDecided: # all agents decided
                break
            
        ys = []
        for n in self.getAllNodes():
            ys.append(self.agents[n].y)
            
        #cascadeCounts = cascadeCounts + [0]*(self.numNodes-len(cascadeCounts))
        if (self.DEBUG): print("Kick sizes: " + str(kickSizes) )
        if (self.DEBUG): print("Cascade log: " + str(cascadeLog) )
        return time, self.countOpinion(1), self.countOpinion(-1), np.mean(ys), np.mean(cascadeCounts), np.std(cascadeCounts), np.mean([abs(k) for k in kickSizes]), self.threshold, sumOfIndividualBayesRisks, cascadeLog
    
        
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



