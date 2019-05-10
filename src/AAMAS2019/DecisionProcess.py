'''
Created on 27 Jan 2017

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

import sys
import os
import configparser
import numpy.random as rand
from matplotlib.backends.backend_pdf import PdfPages
from AAMAS2019 import DecisionNet
from AAMAS2019.MyTypes import NetworkType, DecisionModel, UpdateModel

DEBUG=True

DEFAULT_PROPERTIES_FILENAME = "/Users/joefresna/DecisionsOnNetworks/conf/AAMAS2019.config"

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
    outputTxtFile = config.get('experiment', 'outputTxtFile')
    outputPdfFile = config.get('experiment', 'outputPdfFile')
    cluster = config.getboolean('experiment', 'cluster')
    if (cluster): DEBUG=False
    ## -- Agent params
    decModelStr = config.get('Agent', 'decisionModel')
    if (decModelStr == 'conf-perfect'):
        decModel = DecisionModel.CONFIDENCE
    elif  (decModelStr == 'majority-rand'):
        decModel = DecisionModel.MAJORITY_RAND
    elif  (decModelStr == 'belief-init'): 
        decModel = DecisionModel.BELIEF
    else:
        print("Non valid input for parameter [Agent].decisionModel. Valid values are: 'conf-perfect', 'majority-rand', 'belief-init' ")
        sys.exit()
    maxLoops = config.getint('Agent', 'max_iterations')
    updateConfStr = config.get('Agent', 'updateConf')
    if (updateConfStr == 'no-up'):
        updateConf = UpdateModel.NO_UPDATE
    elif  (updateConfStr == 'optim-up'):
        updateConf = UpdateModel.OPTIMAL
    elif  (updateConfStr == 'belief-up'):
        updateConf = UpdateModel.BELIEF_UP
    elif  (updateConfStr == 'finite-time'):
        updateConf = UpdateModel.FINITE_TIME
    else:
        print("Non valid input for parameter [Agent].updateConf. Valid values are: 'no-up', 'optim-up', 'belief-up', 'finite-time' ")
        sys.exit()
    if (updateConf == UpdateModel.BELIEF_UP or updateConf == UpdateModel.FINITE_TIME) :
        beliefEpsilon = config.getfloat('Agent', 'beliefEpsilon')
        if (updateConf == UpdateModel.FINITE_TIME) :
            finiteTimeExponent = config.getfloat('Agent', 'finiteTimeExponent')
    else:
        beliefEpsilon = None
    ## -- SimpleAgent params
    accuracyMean = config.getfloat('SimpleAgent', 'accuracyMean')
    accuracyStdDev = config.getfloat('SimpleAgent', 'accuracyStdDev')
    truncatePoor = config.getboolean('SimpleAgent', 'truncatePoor')
    ## -- Network params
    numOfNodes = config.getint('Network', 'number_of_nodes')
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
    else:
        print("Non valid input for parameter [Network].netType. Error on value '" + netTypeStr + "'. Valid values are: 'full', 'erdos-renyi', 'barabasi-albert', 'space'")
        sys.exit()
    if (netType == NetworkType.ERSOS_RENYI):
        linkProbability  = config.getfloat('Network', 'link_probability')
    elif (netType == NetworkType.BARABASI_ALBERT):
        numOfEdges  = config.getint('Network', 'number_of_edges')
        if (numOfEdges >= numOfNodes):
            print("STOPPING PROCESS! With barabasi-albert networks the number of edges must be smaller than number of edges. Invalid parameterisation.")
            sys.exit()
    elif (netType == NetworkType.SPACE or netType == NetworkType.SOFT_RGG or 
          netType == NetworkType.FROM_FILE or netType == NetworkType.FROM_FILE_FIXED_COMM ):
        areaSize  = config.getfloat('Network', 'area_size')
        periodicBound = config.getboolean('Network', 'periodic')
        if (netType == NetworkType.SPACE or netType == NetworkType.SOFT_RGG or netType == NetworkType.FROM_FILE_FIXED_COMM):
            commRadius  = config.getfloat('Network', 'communication_radius')
        if (netType == NetworkType.FROM_FILE or netType == NetworkType.FROM_FILE_FIXED_COMM):
            coordinatesFile  = config.get('Network', 'coordinates_file')
            
    dynamicNetwork = config.getboolean('Move', 'dynamic')
    if dynamicNetwork:
        agentSpeed = config.getfloat('Move', 'speed')
        changeDirectionDegStddev = config.getfloat('Move', 'change_direction_deg_stddev')
        plotTrajectory = config.getboolean('Move', 'plotTrajectory')
    else:
        plotTrajectory = False
    
    if DEBUG:
        print( "randomSeed: " + str(randomSeed) )
        print( "numberOfExperiments: " + str(numberOfExperiments) )
        print( "outputTxtFile: " + str(outputTxtFile) )
        print( "outputPdfFile: " + str(outputPdfFile) )
        print( "cluster: " + str(cluster) )
        print( "decModel: " + str(decModelStr) )
        print( "updateModel: " + str(updateConfStr) )
        print( "accuracyMean: " + str(accuracyMean) )
        print( "accuracyStdDev: " + str(accuracyStdDev) )
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
    line = 'seed \t exp \t iter \t pos \t neg \t deg \t clust \t dstd \t cstd  \t conf \t confsd ' + extraInfo + '\n'
    outFile.write(line)
    
    for exp in range(1,numberOfExperiments+1):
        seed = randomSeed*exp
        rand.seed(seed)
        
        ## init the DecNet object
        decNet = DecisionNet.DecNet(netType, decModel, updateConf, numOfNodes, seed, DEBUG)
        
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
        decNet.setSimpleAgent(accuracyMean, accuracyStdDev, truncatePoor)
        if (DEBUG): print("Initialising agents")
        agentParams = []
        if (updateConf == UpdateModel.BELIEF_UP or updateConf == UpdateModel.FINITE_TIME) :
            agentParams.append( beliefEpsilon )
            if (updateConf == UpdateModel.FINITE_TIME) :
                agentParams.append( finiteTimeExponent )
        decNet.initAgents(agentParams)
        
        # Init the agents with their opinion
        decNet.initDecisions()
        
        pp = None
        if not cluster:
            pp = PdfPages(outputPdfFile + '-exp' + str(exp) + '.pdf')
     
        ### Now, the nodes interact with each other
        results = decNet.collectiveDecision(maxLoops, plot=pp, plotTrajectory=plotTrajectory )
        if not cluster:
            pp.close()
        
        line = str(seed) + '\t' + str(exp)
        for res in results:
            line += '\t' + str(res)
        line += '\n'
        outFile.write(line)
    
    outFile.close()   
    if DEBUG: print("Process Ended")
    
