'''
Created on 10 Oct 2016

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''

import sys
import os
import configparser
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
import networkx as nx #@UnresolvedImport
import copy
import DriftDiffMod.DDMclass
import DriftDiffMod.Networks
from statsmodels.graphics import plottools
from matplotlib.backends.backend_pdf import PdfPages

DEBUG=False

DEFAULT_PROPERTIES_FILENAME = "/Users/joefresna/DecisionsOnNetworks/conf/DDMNet.config"

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
    ## -- DDM params
    driftDistribution = config.get('DDM', 'driftDistribution')
    baseDrift = config.getfloat('DDM', 'baseDrift')
    if (driftDistribution == 'uniform'):
        randomDriftRangeMin= config.getfloat('DDM', 'randomDriftRangeMin')  # used for uniform distribution
        randomDriftRangeMax= config.getfloat('DDM', 'randomDriftRangeMax')  # used for uniform distribution
    elif (driftDistribution == 'normal'):
        driftStdDev = config.getfloat('DDM', 'driftStdDev') # used for normal distribution
    noiseStdDev = config.getfloat('DDM', 'noiseStdDev')
    DDMstart = 0
    dt = config.getfloat('DDM', 'dt')
    threshold = config.getfloat('DDM', 'threshold')
    ## -- Network params
    numOfNodes = config.getint('Network', 'number_of_nodes')
    numOfEdges  = config.getint('Network', 'number_of_edges')
    linkProbability  = config.getfloat('Network', 'link_probability')
    ## -- Comparison params
    storeBestDDM = config.getboolean('Comparison', 'storeBestDDM')
    storeConfDDM = config.getboolean('Comparison', 'storeConfDDM')
    storeMajorityRand = config.getboolean('Comparison', 'storeMajorityRand')
    storeMajorityBias = config.getboolean('Comparison', 'storeMajorityBias')
    storeMajorityInhib= config.getboolean('Comparison', 'storeMajorityInhib')
    storeCDCI= config.getboolean('Comparison', 'storeCDCI')
    storeFullNetConf = config.getboolean('Comparison', 'storeFullNetConf')
    storeFullNetMaj = config.getboolean('Comparison', 'storeFullNetMaj')
    
    if DEBUG: 
        if (driftDistribution == 'uniform'): print("DDM drifts are in range [" + str(randomDriftRangeMin) + "," + str(randomDriftRangeMax) + "]" )
        elif (driftDistribution == 'normal'): print("DDM drifts are drawn from normal distribution with mean [" + str(baseDrift) + " and std-dev " + str(driftStdDev) + "" )
    
    ## Open output file and create the first line
    os.makedirs(os.path.dirname(outputTxtFile), exist_ok=True)
    if not cluster: os.makedirs(os.path.dirname(outputPdfFile), exist_ok=True)
    outFile = open(outputTxtFile, 'w')
    line = 'seed \t exp \t run \t iter \t pos \t neg'
    if (storeBestDDM):
        line += '\t bestDDM'
    if (storeConfDDM):
        line += '\t confDDM'
    if (storeFullNetConf):
        line += '\t FNconf'
    if (storeFullNetMaj):
        line += '\t FNmaj'
    if (storeMajorityRand):
        line += '\t iterMR \t posMR \t negMR'
    if (storeMajorityBias):
        line += '\t iterMB \t posMB \t negMB'
    if (storeMajorityInhib):
        line += '\t iterMI \t posMI \t negMI' 
    if (storeCDCI):
        line += '\t iterCDCI \t posCDCI \t negCDCI' 
    line += '\n' 
    outFile.write(line)
    
    for exp in range(1,numberOfExperiments+1):
        rand.seed(randomSeed*exp)
        
        ## Generate a network
        g=DriftDiffMod.Networks.powerlaw_cluster(numOfNodes, numOfEdges, 0, linkProbability, randomSeed)
        #nx.draw(g)
        
        ## initialise a list of DDM (one for each network node n)
        ddms = []
        evols = []
        maxDrift = -float('inf')
        for n in g.nodes():
            # randomly select a drift
            if (driftDistribution == 'uniform'):
                driftRate = baseDrift+rand.uniform(randomDriftRangeMin,randomDriftRangeMax)
            elif (driftDistribution == 'normal'):
                driftRate = rand.normal(baseDrift, driftStdDev)
                
            ddms.append(DriftDiffMod.DDMclass.DDM(
                                                  driftRate, 
                                                  noiseStdDev, 
                                                  DDMstart, 
                                                  dt, 
                                                  threshold))
            evols.append([])
            if (ddms[n].drift > maxDrift):
                maxDriftNode = n
                maxDrift = ddms[n].drift
        
        for run in range(1,repetitionsPerDDM+1):
            # Reset all values
            DriftDiffMod.DDMclass.DDM.converged = 0
            DriftDiffMod.DDMclass.DDM.positive = 0
            DriftDiffMod.DDMclass.DDM.negative = 0
            for n in g.nodes():
                ddms[n].reset(DDMstart)
            
            ## Each timestep, loop through all the nodes and updated their state
            t = 0
            while DriftDiffMod.DDMclass.DDM.converged < nx.number_of_nodes(g):
                t += dt
                for n in g.nodes():
                    if (ddms[n].step()):
                        ddms[n].confidence = ddms[n].logOdds(t);
                        #print(t)
                        #print(ddms[n].y)
                    if (not ddms[n].decided):
                        evols[n].append(ddms[n].y)   
            
            ## Plot the DDMs over time
            if not cluster:
                plt.clf()
                plt.axhline(y=1,  xmin=0, xmax=1, ls='dashed', color='r')
                plt.axhline(y=-1, xmin=0, xmax=1, ls='dashed', color='b')
                plt.ylim(-1.1,1.1)
                #colors=plt.cm.rainbow(np.linspace(0,1,nx.number_of_nodes(g)))
                colors=plottools.rainbow(nx.number_of_nodes(g))
                decisionColors=[]
                for n in g.nodes():
                    plt.plot(np.arange(0,len(evols[n])*dt,dt)[0:len(evols[n])], evols[n],color=colors[n],label='Node '+str(n))
                    if (ddms[n].y >= threshold): decisionColors.append('r') 
                    else: decisionColors.append('b')
            
                # Now add the legend with some customizations.
                #plt.legend(loc='upper right', shadow=True)
                axes = plt.gca()
                lgd = plt.legend(bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
            
                pp = PdfPages(outputPdfFile + '-exp' + str(exp) + '-run' + str(run) + '.pdf')
                pp.savefig(bbox_inches='tight')
                
            ### The individual decisions have been made. I store the decision of the best DDM
            if (ddms[maxDriftNode].y > 0):
                bestDDMDecision = 1
            else:
                bestDDMDecision = -1
            
            ### I store the decision of the most confident DDM 
            if (storeConfDDM):
                maxConf = 0
                for n in g.nodes():
                    if (ddms[n].confidence > maxConf):
                        maxConfNode = n
                        maxConf = ddms[n].confidence
                if (ddms[maxConfNode].y > 0):    
                    confidentDDMDecision = 1
                else:
                    confidentDDMDecision = -1

            ### Compute the fully connected network outcomes
            ### Confidence Log-Odds
            if (storeFullNetConf):
                fullNetConfOpinion=0
                for node in nx.nodes(g):
                    fullNetConfOpinion += (ddms[node].y * ddms[node].confidence)
                    #print('FN Node ' + str(node) + ' : ' + str(ddms[node].y) + ' : ' + str(ddms[node].confidence))
                if (fullNetConfOpinion > 0):
                    fullNetConfOpinion = 1
                elif (fullNetConfOpinion < 0):
                    fullNetConfOpinion = -1
            if (storeFullNetMaj):
                fullNetMajOpinion=0
                for node in nx.nodes(g):
                    fullNetMajOpinion += ddms[node].y
                if (fullNetMajOpinion > 0):
                    fullNetMajOpinion = 1
                elif (fullNetMajOpinion < 0):
                    fullNetMajOpinion = -1
        
            ### Store a copy of the iteration-0 network if different integration methods are compared
            if (storeMajorityRand):
                ddmsMR = copy.deepcopy(ddms)
            if (storeMajorityBias):
                ddmsMB = copy.deepcopy(ddms)
            if (storeMajorityInhib):
                ddmsMI = copy.deepcopy(ddms)
            if (storeCDCI):
                ddmsCDCI = copy.deepcopy(ddms)
        
            ################# CONFIDENCE LOG-ODDS #################
            ### Now, the nodes interact with each other 
            count = 0
            maxLoops = 20
            while (DriftDiffMod.DDMclass.DDM.converged > DriftDiffMod.DDMclass.DDM.positive and 
                   DriftDiffMod.DDMclass.DDM.converged > DriftDiffMod.DDMclass.DDM.negative and count < maxLoops):
                count += 1
                
                if not cluster: 
                    ## Draw initial network
                    fig, axs = plt.subplots(1,2)
                    axs[0].set_title('Iteration n.' + str(count))
                    nx.draw(g,nx.circular_layout(g),node_color=decisionColors,with_labels=True,ax=axs[0])
            
                # create a deep-copy of the current network (for synchronous communication)
                tmp_ddms = copy.deepcopy(ddms)
                # each node updates its opinion as the sum of its weighted opinion plus the neighbours' weighted opinions
                for node in nx.nodes(g):
                    if (DEBUG): print('node ' + str(node))
                    ddms[node].integrateInfoFromNeighbours(tmp_ddms, nx.all_neighbors(g,node), DEBUG)
                    if (DEBUG): print('My new "confidence" (node:' + str(node) + ') is ' + str(ddms[node].confidence))
            
                if (DEBUG): print("Positive nodes: " + str(DriftDiffMod.DDMclass.DDM.positive) + " negative nodes: " + str(DriftDiffMod.DDMclass.DDM.negative) )
            
                if not cluster:
                    ## Draw final network
                    # update colors (red is opt+1, blue is opt-1) and plot
                    decisionColors = []
                    for n in g.nodes():
                        if (ddms[n].y >= threshold):
                            decisionColors.append('r') 
                        else:
                            decisionColors.append('b')
                    nx.draw(g,nx.circular_layout(g),node_color=decisionColors,with_labels=True,ax=axs[1])
            
                ## Print out the report
                if (DEBUG):
                    for n in g.nodes():
                        print('Node ' + str(n) + ' : ' + str(ddms[n].y) + ' : ' + str(ddms[n].confidence))
                
                if not cluster:
                    pp.savefig()
                    plt.close(fig)
            
            numPositive = DriftDiffMod.DDMclass.DDM.positive
            numNegative = DriftDiffMod.DDMclass.DDM.negative 
            
            ################# MAJORITY RULE RANDOM #################
            ## Reset variables and repeat the loop for MajorityRand
            if (storeMajorityRand):
                countMR = 0
                DriftDiffMod.DDMclass.DDM.resetGlobalCounters(ddmsMR)
                for ddm in ddmsMR:
                    ddm.confidence = 1
                    ddm.method = DriftDiffMod.DDMclass.Method.MAJORITY_RAND
                while (DriftDiffMod.DDMclass.DDM.converged > DriftDiffMod.DDMclass.DDM.positive and 
                       DriftDiffMod.DDMclass.DDM.converged > DriftDiffMod.DDMclass.DDM.negative and countMR < maxLoops):
                    countMR += 1
                
                    # create a deep-copy of the current network (for synchronous communication)
                    tmp_ddms = copy.deepcopy(ddmsMR)
                    # each node updates its opinion as the sum of its weighted opinion plus the neighbours' weighted opinions
                    for node in nx.nodes(g):
                        if (DEBUG): print('MR - node ' + str(node))
                        ddmsMR[node].integrateInfoFromNeighbours(tmp_ddms, nx.all_neighbors(g,node), DEBUG)
                
                    if (DEBUG): print("MR - Positive nodes: " + str(DriftDiffMod.DDMclass.DDM.positive) + " negative nodes: " + str(DriftDiffMod.DDMclass.DDM.negative) )
                
                    ## Print out the report
                    if (DEBUG):
                        for n in g.nodes():
                            print('MR - Node ' + str(n) + ' : ' + str(ddmsMR[n].y) + ' : ' + str(ddmsMR[n].confidence))
                
                numPositiveMR = DriftDiffMod.DDMclass.DDM.positive
                numNegativeMR = DriftDiffMod.DDMclass.DDM.negative
                
            ################# MAJORITY RULE BIASED #################
            ## Reset variables and repeat the loop for MajorityRand
            if (storeMajorityBias):
                countMB = 0
                DriftDiffMod.DDMclass.DDM.resetGlobalCounters(ddmsMB)
                for ddm in ddmsMB:
                    ddm.confidence = 1
                    ddm.method = DriftDiffMod.DDMclass.Method.MAJORITY_BIAS
                while (DriftDiffMod.DDMclass.DDM.converged > DriftDiffMod.DDMclass.DDM.positive and 
                       DriftDiffMod.DDMclass.DDM.converged > DriftDiffMod.DDMclass.DDM.negative and countMB < maxLoops):
                    countMB += 1
                
                    # create a deep-copy of the current network (for synchronous communication)
                    tmp_ddms = copy.deepcopy(ddmsMB)
                    # each node updates its opinion as the sum of its weighted opinion plus the neighbours' weighted opinions
                    for node in nx.nodes(g):
                        if (DEBUG): print('MB [it-' + str(countMB) + '] - node ' + str(node))
                        ddmsMB[node].integrateInfoFromNeighbours(tmp_ddms, nx.all_neighbors(g,node), DEBUG)
                
                    if (DEBUG): print("MB - Positive nodes: " + str(DriftDiffMod.DDMclass.DDM.positive) + " negative nodes: " + str(DriftDiffMod.DDMclass.DDM.negative) )
                
                    ## Print out the report
                    if (DEBUG):
                        for n in g.nodes():
                            print('MB - Node ' + str(n) + ' : ' + str(ddmsMB[n].y) + ' : ' + str(ddmsMB[n].confidence))
                
                numPositiveMB = DriftDiffMod.DDMclass.DDM.positive
                numNegativeMB = DriftDiffMod.DDMclass.DDM.negative
                
            ################# MAJORITY RULE INHIBITHED #################
            ## Reset variables and repeat the loop for MajorityInhib
            if (storeMajorityInhib):
                countMI = 0
                DriftDiffMod.DDMclass.DDM.resetGlobalCounters(ddmsMI)
                for ddm in ddmsMI:
                    ddm.confidence = 1
                    ddm.method = DriftDiffMod.DDMclass.Method.MAJORITY_INHIB
                while (DriftDiffMod.DDMclass.DDM.converged > DriftDiffMod.DDMclass.DDM.positive and 
                       DriftDiffMod.DDMclass.DDM.converged > DriftDiffMod.DDMclass.DDM.negative and countMI < maxLoops):
                    countMI += 1
                
                    # create a deep-copy of the current network (for synchronous communication)
                    tmp_ddms = copy.deepcopy(ddmsMI)
                    # each node updates its opinion as the sum of its weighted opinion plus the neighbours' weighted opinions
                    for node in nx.nodes(g):
                        if (DEBUG): print('MI [it-' + str(countMI) + '] - node ' + str(node))
                        ddmsMI[node].integrateInfoFromNeighbours(tmp_ddms, nx.all_neighbors(g,node), DEBUG)
                
                    if (DEBUG): print("MI - Positive nodes: " + str(DriftDiffMod.DDMclass.DDM.positive) + " negative nodes: " + str(DriftDiffMod.DDMclass.DDM.negative) )
                
                    ## Print out the report
                    if (DEBUG):
                        for n in g.nodes():
                            print('MI - Node ' + str(n) + ' : ' + str(ddmsMI[n].y) + ' : ' + str(ddmsMI[n].confidence))
                
                numPositiveMI = DriftDiffMod.DDMclass.DDM.positive
                numNegativeMI = DriftDiffMod.DDMclass.DDM.negative 
                
                
            ##################### CDCI  #####################
            ## Reset variables and repeat the loop for CDCI
            if (storeCDCI):
                countCDCI = 0
                DriftDiffMod.DDMclass.DDM.resetGlobalCounters(ddmsCDCI)
                for ddm in ddmsCDCI:
                    ddm.method = DriftDiffMod.DDMclass.Method.CDCI
                while (DriftDiffMod.DDMclass.DDM.converged > DriftDiffMod.DDMclass.DDM.positive and 
                       DriftDiffMod.DDMclass.DDM.converged > DriftDiffMod.DDMclass.DDM.negative and countCDCI < (maxLoops*100)):
                    countCDCI += 1
                
                    # create a deep-copy of the current network (for synchronous communication)
                    tmp_ddms = copy.deepcopy(ddmsCDCI)
                    # each node updates its opinion through interaction with its neighbours
                    for node in nx.nodes(g):
                        if (DEBUG): print('CDCI [it-' + str(countCDCI) + '] - node ' + str(node))
                        ddmsCDCI[node].integrateInfoFromNeighbours(tmp_ddms, nx.all_neighbors(g,node), DEBUG)
                
                    if (DEBUG): print("CDCI - Positive nodes: " + str(DriftDiffMod.DDMclass.DDM.positive) + " negative nodes: " + str(DriftDiffMod.DDMclass.DDM.negative) )
                
                    ## Print out the report
                    if (DEBUG):
                        for n in g.nodes():
                            print('CDCI - Node ' + str(n) + ' : ' + str(ddmsCDCI[n].y) + ' : ' + str(ddmsCDCI[n].confidence))
                
                numPositiveCDCI = DriftDiffMod.DDMclass.DDM.positive
                numNegativeCDCI = DriftDiffMod.DDMclass.DDM.negative 
            
            if not cluster:
                #plt.show()
                pp.close()
            
            line = str(randomSeed*exp) + '\t' + str(exp) + '\t' + str(run) + '\t' + str(count) + '\t' + str(numPositive) + '\t' + str(numNegative)  
            if (storeBestDDM):
                line += '\t' + str(bestDDMDecision) 
            if (storeConfDDM):
                line += '\t' + str(confidentDDMDecision) 
            if (storeFullNetConf):
                line += '\t' + str(fullNetConfOpinion) 
            if (storeFullNetMaj):
                line += '\t' + str(fullNetMajOpinion) 
            if (storeMajorityRand):
                line += '\t' + str(countMR) + '\t' + str(numPositiveMR) + '\t' + str(numNegativeMR)
            if (storeMajorityBias):
                line += '\t' + str(countMB) + '\t' + str(numPositiveMB) + '\t' + str(numNegativeMB)  
            if (storeMajorityInhib):
                line += '\t' + str(countMI) + '\t' + str(numPositiveMI) + '\t' + str(numNegativeMI)  
            if (storeCDCI):
                line += '\t' + str(countCDCI) + '\t' + str(numPositiveCDCI) + '\t' + str(numNegativeCDCI)  
            line += '\n'
            outFile.write(line)
    
    outFile.close()    
    if DEBUG: print("Process Ended")
    
