'''
Created on 27 Jan 2017

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

import sys
import os
import configparser
# import numpy as np
import numpy.random as rand
# import matplotlib.pyplot as plt
# from statsmodels.graphics import plottools
from matplotlib.backends.backend_pdf import PdfPages
from DecNet import DecisionNet
from DecNet.MyTypes import AgentType, DriftDistribution, NetworkType, DecisionModel

DEBUG=False

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
    elif  (decModelStr == 'log-odds'):
        decModel = DecisionModel.LOGODDS
    else:
        print("Non valid input for parameter [Agent].agentType. Valid values are: 'conf-perfect', 'majority-rand', 'majority-bias', 'majority-inhibit', 'best-acc', 'best-conf', 'log-odds'")
        sys.exit()
    ## -- SimpleAgent params
    if (agentType == AgentType.SIMPLE) :
        accuracyMean = config.getfloat('SimpleAgent', 'accuracyMean')
        accuracyStdDev = config.getfloat('SimpleAgent', 'accuracyStdDev')
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
            randomDriftRangeMin= config.getfloat('DDM', 'randomDriftRangeMin')  # used for uniform distribution
            randomDriftRangeMax= config.getfloat('DDM', 'randomDriftRangeMax')  # used for uniform distribution
        elif (driftDistribution == DriftDistribution.NORMAL):
            baseDrift = config.getfloat('DDM', 'baseDrift')
            driftStdDev = config.getfloat('DDM', 'driftStdDev') # used for normal distribution
        elif (driftDistribution == DriftDistribution.FROM_ACCURACY):
            accuracyMean = config.getfloat('SimpleAgent', 'accuracyMean')
            accuracyStdDev = config.getfloat('SimpleAgent', 'accuracyStdDev')
            
        noiseStdDev = config.getfloat('DDM', 'noiseStdDev')
        DDMstart = 0
        dt = config.getfloat('DDM', 'dt')
        threshold = config.getfloat('DDM', 'threshold')
        interrogationTime = config.getfloat('DDM', 'interrogationTime')
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
    else:
        print("Non valid input for parameter [Network].netType. Error on value '" + netTypeStr + "'. Valid values are: 'full', 'erdos-renyi', 'barabasi-albert'")
        sys.exit()
#     if NetworkType.ERSOS_RENYI in netTypes:
    if (netType == NetworkType.ERSOS_RENYI):
        linkProbability  = config.getfloat('Network', 'link_probability')
    elif (netType == NetworkType.BARABASI_ALBERT):
        numOfEdges  = config.getint('Network', 'number_of_edges')
                
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
            print( "interrogationTime: " + str(interrogationTime) )
        print( "numOfNodes: " + str(numOfNodes) )
        print( "netType: " + str(netTypeStr) )
        if (netType == NetworkType.ERSOS_RENYI):
            print( "linkProbability: " + str(linkProbability) )
        elif (netType == NetworkType.BARABASI_ALBERT):
            print( "numOfEdges: " + str(numOfEdges) )
    
    ## Open output file and create the first line
    os.makedirs(os.path.dirname(outputTxtFile), exist_ok=True)
    if not cluster: os.makedirs(os.path.dirname(outputPdfFile), exist_ok=True)
    outFile = open(outputTxtFile, 'w')
    line = 'seed \t exp \t run \t iter \t pos \t neg \n'
    outFile.write(line)
    
    for exp in range(1,numberOfExperiments+1):
        seed = randomSeed*exp
        rand.seed(seed)
        
        ## init the DecNet object
        decNet = DecisionNet.DecNet(netType, agentType, decModel, numOfNodes, seed, DEBUG)
        
        ## Generate a network
        if (decNet.netType == DecisionNet.NetworkType.FULLY_CONNECTED):
            args = []
        elif (decNet.netType == DecisionNet.NetworkType.ERSOS_RENYI):
            args = [linkProbability]
        elif (decNet.netType == DecisionNet.NetworkType.BARABASI_ALBERT):
            args = [numOfEdges]
        decNet.setNetworkParams(args)
        decNet.initNetwork()
        
        ## initialise the agents (one for each network node n)
        if (decNet.agentType == AgentType.SIMPLE ):
            decNet.setSimpleAgent(accuracyMean, accuracyStdDev)
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
            decNet.setDDMAgent(driftDistribution, baseDrift, noiseStdDev, interrogationTime, args)
        if (DEBUG): print("Initialising agents")
        decNet.initAgents()
        
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
         
            ### Now, the nodes interact with each other
            maxLoops = 20
            results = decNet.collectiveDecision(maxLoops, plot=pp )
            if not cluster:
                #plt.show()
                pp.close()
            
            line = str(seed) + '\t' + str(exp) + '\t' + str(run) + '\t' + str(results[0]) + '\t' + str(results[1]) + '\t' + str(results[2]) + '\n'
            outFile.write(line)
    
    outFile.close()    
    if DEBUG: print("Process Ended")
    
