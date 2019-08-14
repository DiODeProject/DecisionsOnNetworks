'''
Created on 27 Jan 2017

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

import sys
import os
import configparser
import json
import numpy as np
import numpy.random as rand
import scipy.optimize
import scipy.integrate as integrate
from scipy.stats import norm
import math
# import matplotlib.pyplot as plt
# from statsmodels.graphics import plottools
from matplotlib.backends.backend_pdf import PdfPages
from DecNet import DecisionNet
from DecNet.MyTypes import AgentType, DriftDistribution, NetworkType, DecisionModel, UpdateModel

DEBUG=True

DEFAULT_PROPERTIES_FILENAME = "/Users/joefresna/DecisionsOnNetworks/conf/DecNet.config"

if __name__ == '__main__':
    if DEBUG: 
        print("Process Started")
    if len(sys.argv)>1 :
        configFile = sys.argv[1]
    else:
        configFile = DEFAULT_PROPERTIES_FILENAME
    
    if not (os.path.isfile(configFile)):
        print("Impossible to open config file: ", configFile)
        sys.exit()
    
    if DEBUG:
        print("Reading properties file: " + configFile)
    # opening the config file
    config = configparser.ConfigParser()
    config.read(configFile)
   
    ### Load parameters from config file
    ## -- Experiment params
    randomSeed = config.getint('experiment', 'randomSeed')
    numberOfExperiments = config.getint('experiment', 'numberOfExperiments')
    repetitionsPerDDM = config.getint('experiment', 'repetitionsPerDDM')
    outputTxtFile = config.get('experiment', 'outputTxtFile')
    outputPdfFile = config.get('experiment', 'outputPdfFile')
    cluster = config.getboolean('experiment', 'cluster')
    if (cluster): DEBUG=False
    ## -- Agent params
    agentTypeStr = config.get('Agent', 'agentType')
    if (agentTypeStr == 'simple'):
        agentType = AgentType.SIMPLE
    elif  (agentTypeStr == 'DDM'):
        agentType = AgentType.DDM
    else:
        print("Non valid input for parameter [Agent].agentType. Valid values are: 'simple', 'DDM'")
        sys.exit()
    decModelStr = config.get('Agent', 'decisionModel')
    if (decModelStr == 'conf-perfect'):
        decModel = DecisionModel.CONFIDENCE
    elif  (decModelStr == 'majority-rand'):
        decModel = DecisionModel.MAJORITY_RAND
    elif  (decModelStr == 'majority-bias'):
        decModel = DecisionModel.MAJORITY_BIAS
    elif  (decModelStr == 'majority-inhibit'):
        decModel = DecisionModel.MAJORITY_INHIB
    elif  (decModelStr == 'best-acc'):
        decModel = DecisionModel.BEST_ACC
    elif  (decModelStr == 'best-conf'):
        decModel = DecisionModel.BEST_CONF
    elif  (decModelStr == 'log-odds-perfect'):
        decModel = DecisionModel.LOGODDS_PERFECT
    elif  (decModelStr == 'log-odds-combo'):
        decModel = DecisionModel.LOGODDS_COMBO
    elif  (decModelStr == 'log-odds-distr'):
        decModel = DecisionModel.LOGODDS_DISTRIBUTION
    elif  (decModelStr == 'log-odds-approx'): 
        decModel = DecisionModel.LOGODDS_APPROX
    elif  (decModelStr == 'belief-init'): 
        decModel = DecisionModel.BELIEF
    else:
        print("Non valid input for parameter [Agent].decisionModel. Valid values are: 'conf-perfect', 'majority-rand', 'majority-bias', 'majority-inhibit', 'best-acc', 'best-conf', 'log-odds-perfect', 'log-odds-combo', 'log-odds-distr', 'log-odds-approx', 'belief-init' ")
        sys.exit()
    maxLoops = config.getint('Agent', 'max_iterations')
    updateConfStr = config.get('Agent', 'updateConf')
    if (updateConfStr == 'no-up'):
        updateConf = UpdateModel.NO_UPDATE
    elif  (updateConfStr == 'theta-up'):
        updateConf = UpdateModel.THETA_UPDATE
    elif  (updateConfStr == 'theta-norm'):
        updateConf = UpdateModel.THETA_NORM
    elif  (updateConfStr == 'optim-up'):
        updateConf = UpdateModel.OPTIMAL
    elif  (updateConfStr == 'belief-up'):
        updateConf = UpdateModel.BELIEF_UP
    elif  (updateConfStr == 'finite-time'):
        updateConf = UpdateModel.FINITE_TIME
    else:
        print("Non valid input for parameter [Agent].updateConf. Valid values are: 'no-up', 'theta-up', 'theta-norm', 'optim-up', 'belief-up', 'finite-time' ")
        sys.exit()
    if (updateConf == UpdateModel.BELIEF_UP or updateConf == UpdateModel.FINITE_TIME) :
        beliefEpsilon = config.getfloat('Agent', 'beliefEpsilon')
        if (updateConf == UpdateModel.FINITE_TIME) :
            finiteTimeExponent = config.getfloat('Agent', 'finiteTimeExponent')
    else:
        beliefEpsilon = None
    ## -- SimpleAgent params
    if (agentType == AgentType.SIMPLE) :
        accuracyMean = config.getfloat('SimpleAgent', 'accuracyMean')
        accuracyStdDev = config.getfloat('SimpleAgent', 'accuracyStdDev')
        truncatePoor = config.getboolean('SimpleAgent', 'truncatePoor')
    ## -- DDM params
    if (agentType == AgentType.DDM) :
        driftDistributionStr = config.get('DDM', 'driftDistribution')
        if (driftDistributionStr == 'uniform'):
            driftDistribution = DriftDistribution.UNIFORM
        elif  (driftDistributionStr == 'normal'):
            driftDistribution = DriftDistribution.NORMAL
        elif  (driftDistributionStr == 'from_accuracy'):
            driftDistribution = DriftDistribution.FROM_ACCURACY
        else:
            print("Non valid input for parameter [DDM].driftDistribution. Valid values are: 'uniform', 'normal', 'from_accuracy'")
            sys.exit()
        if (driftDistribution == DriftDistribution.UNIFORM):
            baseDrift = config.getfloat('DDM', 'baseDrift')
            randomDriftRangePlusMinus= config.getfloat('DDM', 'randomDriftRangePlusMinus')  # used for uniform distribution
            randomDriftRangeMin= baseDrift - abs( randomDriftRangePlusMinus )
            randomDriftRangeMax= baseDrift + abs( randomDriftRangePlusMinus )
        elif (driftDistribution == DriftDistribution.NORMAL):
            baseDrift = config.getfloat('DDM', 'baseDrift')
            driftStdDev = config.getfloat('DDM', 'driftStdDev') # used for normal distribution
        elif (driftDistribution == DriftDistribution.FROM_ACCURACY):
            accuracyMean = config.getfloat('SimpleAgent', 'accuracyMean')
            accuracyStdDev = config.getfloat('SimpleAgent', 'accuracyStdDev')
            
        noiseStdDev = config.getfloat('DDM', 'noiseStdDev')
        dt = config.getfloat('DDM', 'dt')
        threshold = config.getfloat('DDM', 'threshold')
        interrogationTime = config.getfloat('DDM', 'interrogationTime')
        interrogationTimeFromAccuracy = config.getboolean('DDM', 'interrogationTimeFromAccuracy')
        prior = config.getfloat('DDM', 'prior')
        costMatrix = json.loads( config['DDM']['costMatrix'] )
    ## -- Network params
    numOfNodes = config.getint('Network', 'number_of_nodes')
#     netTypesStr = config.get('Network', 'netType').split(',')
#     netTypes = []
#     for netTypeStr in netTypesStr :
#         if (netTypeStr.lstrip() == 'power-law'):
#             netTypes.append(NetworkType.POWER_LAW)
#         elif (netTypeStr.lstrip() == 'erdos-renyi'):
#             netTypes.append(NetworkType.ERSOS_RENYI)
#         elif (netTypeStr.lstrip() == 'full'):
#             netTypes.append(NetworkType.FULLY_CONNECTED)
#         else:
#             print("Non valid input for parameter [Network].netType. Error on value '" + netTypeStr.lstrip() + "'. Valid values are: 'full', 'erdos-renyi', 'power-law'")
#             sys.exit()
    netTypeStr = config.get('Network', 'netType')
    if (netTypeStr == 'full'):
        netType =NetworkType.FULLY_CONNECTED
    elif (netTypeStr == 'erdos-renyi'):
        netType = NetworkType.ERSOS_RENYI
    elif (netTypeStr == 'barabasi-albert'):
        netType = NetworkType.BARABASI_ALBERT
    elif (netTypeStr == 'space'):
        netType = NetworkType.SPACE
    elif (netTypeStr == 'soft-rgg'):
        netType = NetworkType.SOFT_RGG
    elif (netTypeStr == 'from-file'):
        netType = NetworkType.FROM_FILE
    elif (netTypeStr == 'from-file-fixComm'):
        netType = NetworkType.FROM_FILE_FIXED_COMM
    elif (netTypeStr == 'rgg-fixed-degree'):
        netType = NetworkType.RGG_FIXED_DEGREE
    else:
        print("Non valid input for parameter [Network].netType. Error on value '" + netTypeStr + "'. Valid values are: 'full', 'erdos-renyi', 'barabasi-albert', 'space'")
        sys.exit()
#     if NetworkType.ERSOS_RENYI in netTypes:
    if (netType == NetworkType.ERSOS_RENYI):
        linkProbability  = config.getfloat('Network', 'link_probability')
    elif (netType == NetworkType.BARABASI_ALBERT):
        numOfEdges  = config.getint('Network', 'number_of_edges')
        if (numOfEdges >= numOfNodes):
            print("STOPPING PROCESS! With barabasi-albert networks the number of edges must be smaller than number of edges. Invalid parameterisation.")
            sys.exit()
    elif (netType == NetworkType.SPACE or netType == NetworkType.SOFT_RGG or netType == NetworkType.FROM_FILE or 
          netType == NetworkType.FROM_FILE_FIXED_COMM or netType == NetworkType.RGG_FIXED_DEGREE ):
        areaSize  = config.getfloat('Network', 'area_size')
        periodicBound = config.getboolean('Network', 'periodic')
        if (netType == NetworkType.SPACE or netType == NetworkType.SOFT_RGG or netType == NetworkType.FROM_FILE_FIXED_COMM or netType == NetworkType.RGG_FIXED_DEGREE):
            commRadius  = config.getfloat('Network', 'communication_radius')
        if (netType == NetworkType.RGG_FIXED_DEGREE): # convert the fixed degree (now called commRadius) in the actual commRadius
            commRadius = math.sqrt( commRadius * areaSize * areaSize / (math.pi * numOfNodes) )
            netType = NetworkType.SPACE # then it can be treated as normal SPACE networks
        if (netType == NetworkType.FROM_FILE or netType == NetworkType.FROM_FILE_FIXED_COMM):
            coordinatesFile  = config.get('Network', 'coordinates_file')
            
    dynamicNetwork = config.getboolean('Move', 'dynamic')
    if dynamicNetwork:
        agentSpeed = config.getfloat('Move', 'speed')
        changeDirectionDegStddev = config.getfloat('Move', 'change_direction_deg_stddev')
        plotTrajectory = config.getboolean('Move', 'plotTrajectory')
    else:
        plotTrajectory = False
    ## -- Comparison params
#     storeBestDDM = config.getboolean('Comparison', 'storeBestDDM')
#     storeConfDDM = config.getboolean('Comparison', 'storeConfDDM')
#     storeMajorityRand = config.getboolean('Comparison', 'storeMajorityRand')
#     storeMajorityBias = config.getboolean('Comparison', 'storeMajorityBias')
#     storeMajorityInhib= config.getboolean('Comparison', 'storeMajorityInhib')
#     storeCDCI= config.getboolean('Comparison', 'storeCDCI')
#     storeFullNetConf = config.getboolean('Comparison', 'storeFullNetConf')
#     storeFullNetMaj = config.getboolean('Comparison', 'storeFullNetMaj')
    
    if DEBUG:
        print( "randomSeed: " + str(randomSeed) )
        print( "numberOfExperiments: " + str(numberOfExperiments) )
        print( "repetitionsPerDDM: " + str(repetitionsPerDDM) )
        print( "outputTxtFile: " + str(outputTxtFile) )
        print( "outputPdfFile: " + str(outputPdfFile) )
        print( "cluster: " + str(cluster) )
        print( "agentType: " + str(agentTypeStr) )
        print( "decModel: " + str(decModelStr) )
        print( "updateModel: " + str(updateConfStr) )
        if (agentType == AgentType.SIMPLE) :
            print( "accuracyMean: " + str(accuracyMean) )
            print( "accuracyStdDev: " + str(accuracyStdDev) )
        if (agentType == AgentType.DDM) :
            print( "driftDistribution: " + str(driftDistributionStr) )
            if (driftDistribution == DriftDistribution.UNIFORM): print("DDM drifts are in range [" + str(randomDriftRangeMin) + "," + str(randomDriftRangeMax) + "]" )
            elif (driftDistribution == DriftDistribution.NORMAL): print("DDM drifts are drawn from normal distribution with mean " + str(baseDrift) + " and std-dev " + str(driftStdDev) + "" )
            print( "noiseStdDev: " + str(noiseStdDev) )
            print( "dt: " + str(dt) )
            print( "threshold: " + str(threshold) )
            print( "interrogationTimeFromAccuracy: " + str(interrogationTimeFromAccuracy) )
            print( "interrogationTime: " + str(interrogationTime) )
            print( "prior: " + str(prior) )
            print( "costMatrix: " + str(costMatrix) )
        print( "numOfNodes: " + str(numOfNodes) )
        print( "netType: " + str(netTypeStr) )
        if (netType == NetworkType.ERSOS_RENYI):
            print( "linkProbability: " + str(linkProbability) )
        elif (netType == NetworkType.BARABASI_ALBERT):
            print( "numOfEdges: " + str(numOfEdges) )
        elif (netType == NetworkType.SPACE or netType == NetworkType.SOFT_RGG):
            print( "areaSize: " + str(areaSize) )
            print( "commRadius: " + str(commRadius) )
            print( "Periodic boundary conditions: " + str(periodicBound) )
        elif (netType == NetworkType.FROM_FILE):
            print( "areaSize: " + str(areaSize) )
        elif (netType == NetworkType.FROM_FILE_FIXED_COMM):
            print( "areaSize: " + str(areaSize) )
            print( "commRadius: " + str(commRadius) )
        print( "Dynamic network: " + str(dynamicNetwork) )
        if dynamicNetwork:
            print( "Agent speed: " + str(agentSpeed) )
            print( "Agent change of direction std-dev: " + str(changeDirectionDegStddev) + " degrees")
            print( "Plot trajectory: " + str(plotTrajectory) )
    
    ## Open output file and create the first line
    os.makedirs(os.path.dirname(outputTxtFile), exist_ok=True)
    if not cluster: os.makedirs(os.path.dirname(outputPdfFile), exist_ok=True)
    outFile = open(outputTxtFile, 'w')
    extraInfo = ''
#     if (netType == DecisionNet.NetworkType.FULLY_CONNECTED):
#         extraInfo = '\t acc ' # for fully-connected network we also estimate the expected group accuracy
    line = 'seed \t exp \t run \t iter \t pos \t neg \t deg \t clust \t dstd \t cstd  \t conf \t confsd ' + extraInfo + '\n'
    outFile.write(line)
    
    for exp in range(1,numberOfExperiments+1):
        seed = randomSeed*exp
        rand.seed(seed)
        
        ## init the DecNet object
        decNet = DecisionNet.DecNet(netType, agentType, decModel, updateConf, numOfNodes, seed, DEBUG)
        decNet.logFilename = "/Users/joefresna/DecisionsOnNetworks/data/conf-log" + str(exp) + ".txt"
        
        ## Generate a network
        if (decNet.netType == NetworkType.FULLY_CONNECTED):
            args = []
        elif (decNet.netType == NetworkType.ERSOS_RENYI):
            args = [linkProbability]
        elif (decNet.netType == NetworkType.BARABASI_ALBERT):
            args = [numOfEdges]
        elif (decNet.netType == NetworkType.SPACE or netType == NetworkType.SOFT_RGG):
            args = [areaSize, commRadius, periodicBound]
        elif (decNet.netType == NetworkType.FROM_FILE):
            args = [areaSize, coordinatesFile]
        elif (decNet.netType == NetworkType.FROM_FILE_FIXED_COMM):
            args = [areaSize, coordinatesFile, commRadius]
        if dynamicNetwork:
            args.extend([agentSpeed, changeDirectionDegStddev])
        decNet.setNetworkParams(dynamicNetwork, args)
        decNet.initNetwork()
        
        ## initialise the agents (one for each network node n)
        if (decNet.agentType == AgentType.SIMPLE ):
            decNet.setSimpleAgent(accuracyMean, accuracyStdDev, truncatePoor)
        elif (decNet.agentType == AgentType.DDM ):
            args = []
            if (driftDistribution == DriftDistribution.UNIFORM):
                args.append( randomDriftRangeMin ) 
                args.append( randomDriftRangeMax )
            elif (driftDistribution == DriftDistribution.NORMAL):
                args.append( driftStdDev )
            elif (driftDistribution == DriftDistribution.FROM_ACCURACY):
                baseDrift = 0
                args.append( accuracyMean )
                args.append( accuracyStdDev )
            if interrogationTimeFromAccuracy:
                args.append(costMatrix)
                def integrationFuncDdmThreshNormal(drift, costMatrix, noiseStdDev, baseDrift, driftStdDev):
                    prob = norm.pdf(drift, baseDrift, driftStdDev)
                    drift = abs(drift)
                    def compute_BR_opt_thresh(thresh): # from Eq. (5.6) of Bogacz et al. Psy.Rev. 2016
                        return costMatrix[1]/costMatrix[0] * 2 * (drift**2) / (noiseStdDev**2) - 4 * drift * thresh / (noiseStdDev**2) + np.exp( - 2 * drift * thresh / (noiseStdDev**2) ) -  np.exp( 2 * drift * thresh / (noiseStdDev**2) )             
                    pdf = scipy.optimize.fsolve( compute_BR_opt_thresh, 0.25*drift*costMatrix[1]/costMatrix[0], maxfev=1000 )[0] # second parameter is the starting point which is set to the approximation of Eq. (5.7) of Bogacz et al. Psy.Rev. 2016
                    return pdf*prob
                integrResult = integrate.quad( integrationFuncDdmThreshNormal, -np.inf, np.inf, args=( costMatrix, noiseStdDev, baseDrift, driftStdDev ) )
                thresh = integrResult[0]
                expectedAccuracy = 1 - ( 1.0/ ( 1 + np.exp(2*baseDrift*thresh/(noiseStdDev**2)) ))
                interrogationTime = expectedAccuracy
            decNet.setDDMAgent(driftDistribution, baseDrift, noiseStdDev, interrogationTimeFromAccuracy, interrogationTime, prior, args)
        if (DEBUG): print("Initialising agents")
        agentParams = []
        if (updateConf == UpdateModel.BELIEF_UP or updateConf == UpdateModel.FINITE_TIME) :
            agentParams.append( beliefEpsilon )
            if (updateConf == UpdateModel.FINITE_TIME) :
                agentParams.append( finiteTimeExponent )
        decNet.initAgents(agentParams)
        
        for run in range(1,repetitionsPerDDM+1):
            # Init the agents with their opinion
            decNet.initDecisions()
            
#             ## Plot the DDMs over time
#             if not cluster:
#                 plt.clf()
#                 plt.axhline(y=1,  xmin=0, xmax=1, ls='dashed', color='r')
#                 plt.axhline(y=-1, xmin=0, xmax=1, ls='dashed', color='b')
#                 plt.ylim(-1.1,1.1)
#                 colors=plottools.rainbow(len(decNet.getAllNodes()))
#                 decisionColors=[]
#                 for n in decNet.getAllNodes():
#                     plt.plot(np.arange(0,len(evols[n])*dt,dt)[0:len(evols[n])], evols[n],color=colors[n],label='Node '+str(n))
#                     if (decNet.agents[n].opinion > 0): decisionColors.append('r') 
#                     elif  (decNet.agents[n].opinion < 0): decisionColors.append('b')
#                     else: decisionColors.append('w')
#              
#                 # Now add the legend with some customizations.
#                 plt.legend(loc='upper right', shadow=True)
#                 axes = plt.gca()
#                 lgd = plt.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
            pp = None
            if not cluster:
                pp = PdfPages(outputPdfFile + '-exp' + str(exp) + '-run' + str(run) + '.pdf')
#                 pp.savefig(bbox_inches='tight')
                #trajectoryPlotPrefix = outputPdfFile + '-exp' + str(exp) + '-run' + str(run) if plotTrajectory else None 
         
            ### Now, the nodes interact with each other
            results = decNet.collectiveDecision(maxLoops, plot=pp, plotTrajectory=plotTrajectory )
            if not cluster:
                #plt.show()
                pp.close()
            
#             extraInfo = ''
#             if (decNet.netType == DecisionNet.NetworkType.FULLY_CONNECTED):
#                 extraInfo = '\t' + str(results[3]) # for fully-connected network we also estimate the expected group accuracy
#             line = str(seed) + '\t' + str(exp) + '\t' + str(run) + '\t' + str(results[0]) + '\t' + str(results[1]) + '\t' + str(results[2]) + extraInfo + '\n'
            #line = str(seed) + '\t' + str(exp) + '\t' + str(run) + '\t' + str(results[0]) + '\t' + str(results[1]) + '\t' + str(results[2]) + '\t' + str(results[3]) + '\t' + str(results[4]) + '\t' + str(results[5]) + '\t' + str(results[6]) + '\n'
            line = str(seed) + '\t' + str(exp) + '\t' + str(run)
            for res in results:
                line += '\t' + str(res)
            line += '\n'
            outFile.write(line)
    
    outFile.close()   
    if DEBUG: print("Process Ended")
    
    
