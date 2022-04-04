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
import scipy.integrate as integrate
from scipy.stats import norm
import scipy.optimize
import warnings

def computeDriftFromAccuracy(accuracy, noise, interrogationTime):
    return math.sqrt(2) * noise * scipy.special.erfcinv(2 - 2*accuracy)/math.sqrt(interrogationTime)

def infiniteSeriesBorodinSalminen(t, u, v, K):
    summedValue = 0
    for k in K:
        summedValue += (v -u +2*k*v) * np.exp( -(v -u +2*k*v)**2 / (2*t) ) / np.sqrt(2*np.pi*(t**3))
    return summedValue

def integrationFuncDriftNormalFreeResp(drift, thresh, noise, time, sign, popMean, popStdDev, misinformed, normalisationOfNormProbForResampling):
    if misinformed:
        pdf_normal_dist = norm.pdf(drift, sign*popMean, popStdDev)
        signedDrift = drift
    else:
        pdf_normal_dist = norm.pdf(drift, popMean, popStdDev) / normalisationOfNormProbForResampling
        signedDrift = sign*abs(drift)
    expected_accuracy = 1- (1 / (1 + np.exp( 2*signedDrift*thresh / noise**2 ) ))
         
    ## to compute the pdf_first_passage_time for DDM with two thresholds we use Eq. (14) from https://arxiv.org/pdf/1508.03373.pdf
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        try:
            denominator_pdf_first_passage_time = 1 - (1 - np.exp(-2*signedDrift*thresh/(noise**2)) ) / ( np.exp(2*signedDrift*thresh/(noise**2)) - np.exp(-2*signedDrift*thresh/(noise**2)) )
            #if sign==-1: print("The denominator is " + str(denominator_pdf_first_passage_time))
            if math.isnan(denominator_pdf_first_passage_time ): print(signedDrift)
            pdf_first_passage_time = np.exp( -(signedDrift**2)*time/(2*noise**2) + signedDrift*thresh/(noise**2) ) * infiniteSeriesBorodinSalminen(time, thresh/noise, 2*thresh/noise, np.arange(-10,11)) / denominator_pdf_first_passage_time
        except RuntimeWarning:
            pdf_first_passage_time = 0
    #print("The FPT is " + str(pdf_first_passage_time))
    
    ## equation for a single threshold drift
    #signedDrift = signedDrift/noise
    #thresh = thresh/noise
    #pdf_first_passage_time = (thresh / np.sqrt( 2 * np.pi * time**3) ) * np.exp( (-(thresh - signedDrift * time)**2) / (2 * time))
              
    return pdf_normal_dist* expected_accuracy * pdf_first_passage_time

def estimateConfFromDistribution(signedThreshold, time, noise, prior, meanDrift, stdDrift, misinformed):
    integrationLowerLimit = -np.inf if misinformed else 0
    normalisationOfNormProbForResampling = 1 if misinformed else integrate.quad( norm.pdf, 0, np.inf, args=(meanDrift, stdDrift))[0]
    numerator   = integrate.quad( integrationFuncDriftNormalFreeResp, integrationLowerLimit, np.inf, args=( abs(signedThreshold), noise, time, +1, meanDrift, stdDrift, misinformed, normalisationOfNormProbForResampling ) )
    denominator = integrate.quad( integrationFuncDriftNormalFreeResp, integrationLowerLimit, np.inf, args=( abs(signedThreshold), noise, time, -1, meanDrift, stdDrift, misinformed, normalisationOfNormProbForResampling ) )
    priorComponent = 0 if (prior == 0.5) else math.log( prior/ (1-prior))  # if the sing of signedThreshold is negative, through math.copysign I change the sing of the priorComponent (which is equivalent to power to -1 the log argument)
    return math.log( numerator[0]/denominator[0] ) + math.copysign(priorComponent, signedThreshold)

def integrFunc(u):
    val = np.exp( -(u**2)/2) / math.sqrt(2 * math.pi)
    return val 

def accuracyInterrogation(drift, time, noise):
    return 1 - integrate.quad( integrFunc, -np.inf, -(drift * math.sqrt(time))/noise )[0]

def integrationFuncAccuracyInterrogation(drift, time, noise, popMean, popStdDev, misinformed):        
    pdf_normal_dist = norm.pdf(drift, popMean, popStdDev)
    if not misinformed: drift = abs(drift)
    expected_accuracy = accuracyInterrogation(drift, time, noise) 
    return pdf_normal_dist* expected_accuracy

def accuracyFreeResponse(drift, thresh, noise):
    return 1- (1 / (1 + np.exp( 2*drift*thresh / (noise**2) ) ))

def integrationFuncAccuracyFreeResponse(drift, thresh, noise, popMean, popStdDev, misinformed):
    pdf_normal_dist = norm.pdf(drift, popMean, popStdDev)
    if not misinformed: drift = abs(drift)
    expected_accuracy = accuracyFreeResponse(drift, thresh, noise) 
    return pdf_normal_dist* expected_accuracy

def debugConfidenceKick():
#     noiseStdDev = 1
#     thresh = 1
#     prior = 0.5
#     meanDrift = 0.5
#     stdDrift = 0.1
    noiseStdDev = 0.5
    thresh = 1
    prior = 0.5
    meanDrift = 0.05
    stdDrift = 0.1
    misinformed = False
    
    tmpTime = 1.1 
    tmpDrift=0.5
    print("pdf_normal_dist = " + str(norm.pdf(tmpDrift, meanDrift, stdDrift)))
    print("expected_accuracy = " + str(1- (1 / (1 + np.exp( 2*tmpDrift*thresh / noiseStdDev**2 ) ))))
    print("pdf_first_passage_time = " + str((thresh / np.sqrt( 2 * np.pi * tmpTime**3) ) * np.exp( (-(thresh - tmpDrift * tmpTime)**2) / (2 * tmpTime)))) 
    
    tmp = integrate.quad( integrationFuncDriftNormalFreeResp, -np.inf, np.inf, args=( thresh, noiseStdDev, tmpTime, +1, meanDrift, stdDrift, misinformed ) )
    print("for time " + str(tmpTime) + " pdf of fpt for plus thresh is " + str(tmp[0]))
    tmp = integrate.quad( integrationFuncDriftNormalFreeResp, -np.inf, np.inf, args=( -thresh, noiseStdDev, tmpTime, +1, meanDrift, stdDrift, misinformed ) )
    print("for time " + str(tmpTime) + " pdf of fpt for minus thresh is " + str(tmp[0]))
        
    out = '{'
    #for time in [0.05, 0.1, 0.2, 1, 2, 3]:
    for time in np.arange(0.1, 10, 1):
        conf = estimateConfFromDistribution(thresh, time, noiseStdDev, prior, meanDrift, stdDrift, misinformed)
        out += "{" + str(time) + "," + str(conf) + "},"
    print(out[:-1] + '}')
    print('Done')

def integrationFuncDdmThreshNormal(drift, costMatrix, noiseStdDev, baseDrift, driftStdDev, normalisationOfNormProbForResampling):
    prob = norm.pdf(drift, baseDrift, driftStdDev) / normalisationOfNormProbForResampling
    #if not misinformed: drift = abs(drift)
    def compute_BR_opt_thresh(thresh): # from Eq. (5.6) of Bogacz et al. Psy.Rev. 2016
        return costMatrix[1]/costMatrix[0] * 2 * (drift**2) / (noiseStdDev**2) - 4 * drift * thresh / (noiseStdDev**2) + np.exp( - 2 * drift * thresh / (noiseStdDev**2) ) -  np.exp( 2 * drift * thresh / (noiseStdDev**2) )             
    pdf = scipy.optimize.fsolve( compute_BR_opt_thresh, 0.25*drift*costMatrix[1]/costMatrix[0], maxfev=5000 )[0] # second parameter is the starting point which is set to the approximation of Eq. (5.7) of Bogacz et al. Psy.Rev. 2016
    return pdf*prob

def computedThreshold(costMatrix, noiseStdDev, meanDrift, driftStdDev, misinformed):
    integrationLowerLimit = -np.inf if misinformed else 0
    normalisationOfNormProbForResampling = 1 if misinformed else integrate.quad( norm.pdf, 0, np.inf, args=(meanDrift, driftStdDev))[0]
    threshold = integrate.quad( integrationFuncDdmThreshNormal, integrationLowerLimit, np.inf, args=( costMatrix, noiseStdDev, meanDrift, driftStdDev, normalisationOfNormProbForResampling) )
    #print('thresh is ' + str(threshold))
    return threshold[0]

def debugAccuracyFreeResponse():
    noiseStdDev = 0.5
    meanDrift = 0.05
    for meanDrift in [0.05, 0.1, 0.2]:
        print("mean drift " + str(meanDrift))
        costMatrix = [1, 100]
        misinformed = False
    #     thresh = 0.4122258202119302
    #     acc = accuracyFreeResponse(meanDrift, thresh, noiseStdDev)
    #     print( "Mean accuracy: " + str(acc) )
        out = 'c('
        stdDeviations = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
        for stdDrift in stdDeviations:
            thresh = computedThreshold(costMatrix, noiseStdDev, meanDrift, stdDrift, misinformed)
            integr = integrate.quad( integrationFuncAccuracyFreeResponse, -np.inf, np.inf, args=( thresh, noiseStdDev, meanDrift, stdDrift, misinformed ) )
            #print( "Integrated accuracy: " + str(integr) )
            out += str(integr[0]) + ","
            #out += str(thresh) + ","
        print('devsFR <- c(' + str( stdDeviations )[1:-1] + ')')
        print('valsFR <- ' + out[:-1] + ')')

def computedInterrogationTime(costMatrix, noise, meanDrift, driftStdDev, misinformed):
    integrResult = integrate.quad( integrationFuncDdmThreshNormal, -np.inf, np.inf, args=( costMatrix, noise, meanDrift, driftStdDev, misinformed ) )
    thresh = integrResult[0]
    def integrationFuncExpAccuracy(drift, thresh, noiseStdDev, baseDrift, driftStdDev, misinformed ):
        prob = norm.pdf(drift, baseDrift, driftStdDev)
        if not misinformed: drift = abs(drift)
        expAcc = 1 - ( 1.0/ ( 1 + np.exp(2*drift*thresh/(noiseStdDev**2)) ))
        return prob*expAcc
    #expectedAcc = 1 - ( 1.0/ ( 1 + np.exp(2*meanDrift*thresh/(noise**2)) ))
    expectedAcc = integrate.quad( integrationFuncExpAccuracy, -np.inf, np.inf, args=( thresh, noise, meanDrift, driftStdDev, misinformed ) )[0]
    interrogationTime = costMatrix[1]/costMatrix[0] * (2*expectedAcc - 1) * math.log( expectedAcc / (1-expectedAcc) ) / ( 2*math.log(expectedAcc / (1-expectedAcc)) - 1/expectedAcc + 1/(1-expectedAcc) )
    return interrogationTime

def debugAccuracyInterrogation():
    noiseStdDev = 0.5
    costMatrix = [1, 20]
    misinformed = False
    #meanDrift = 0.05
    #time = 0.2492898885850068   #0.5938547459929719 #0.9607487316511182
    #acc = accuracyInterrogation(meanDrift, time, noiseStdDev)
    #print( "Mean accuracy: " + str(acc) )
    for meanDrift in [0.05, 0.1, 0.2]:
        print("mean drift " + str(meanDrift))
        out = 'c('
        stdDeviations = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
        for stdDrift in stdDeviations:
            time = computedInterrogationTime(costMatrix, noiseStdDev, meanDrift, stdDrift, misinformed)
            integr = integrate.quad( integrationFuncAccuracyInterrogation, -np.inf, np.inf, args=( time, noiseStdDev, meanDrift, stdDrift, misinformed ) )
            #print( "Integrated accuracy: " + str(integr) )
            out += str(integr[0]) + ","
        print('devsIn <- c(' + str( stdDeviations )[1:-1] + ')')
        print('valsIn <- ' + out[:-1] + ')')
    

if __name__ == '__main__':
    costMatrix = [1, 100]
    noiseStdDev = 1
    prior = 0.5
    for meanDrift in [0.2]: #[0.05, 0.1, 0.2]:
        print("mean drift " + str(meanDrift))
        misinformed = False
        out = '{'
        stdDeviations = [0.5,0.75,1]
        #stdDeviations = [0.2]
        for stdDrift in stdDeviations:
            thresh = computedThreshold(costMatrix, noiseStdDev, meanDrift, stdDrift, misinformed)
            xs = []
            ys = []
            for time in np.arange(0.25, 20, 0.25):
                if time==0: time+=0.0000000000000000001
                kickSize = estimateConfFromDistribution(+1, time, noiseStdDev, prior, meanDrift, stdDrift, misinformed)
                out+='{'+str(stdDrift)+","+str('{:.20f}'.format(time))+","+str(kickSize/thresh)+"},"

        out=out[:-1]+'}'
        print(out)
    
if __name__ == '__main2__':
    costMatrix = [1, 100]
    noiseStdDev = 1
    prior = 0.5
    for meanDrift in [0.2]: #[0.05, 0.1, 0.2]:
        print("mean drift " + str(meanDrift))
        misinformed = False
        out = '{'
        out2 = '{'
        stdDeviations = np.arange(0.05, 0.51, 0.05)
        #stdDeviations = [0.2]
        for stdDrift in stdDeviations:
            thresh = computedThreshold(costMatrix, noiseStdDev, meanDrift, stdDrift, misinformed)
            xs = []
            ys = []
            out2 += '{'+str(stdDrift)+","+str(thresh)+'},'
#             minRange = max(0.0001, meanDrift-3*stdDrift) if not misinformed else meanDrift-3*stdDrift
#             stepSize = 6*stdDrift/20
#            maxRange = meanDrift+3*stdDrift+stepSize/2
            stepSize = min(stdDeviations)*6/10
            minRange = max(0,meanDrift-3*max(stdDeviations))
            maxRange = meanDrift+3*max(stdDeviations)+stepSize/2
            for drift in np.arange(minRange, maxRange, stepSize):
                if drift==0: drift+=0.0000000000000000001
                expectedConvTime = (thresh/drift) * np.tanh( drift * thresh/ (noiseStdDev**2) )
                kickSize = estimateConfFromDistribution(+1, expectedConvTime, noiseStdDev, prior, meanDrift, stdDrift, misinformed)
                out+='{'+str(stdDrift)+","+str('{:.20f}'.format(drift))+","+str(kickSize/thresh)+"},"
#                xs.append(drift)
#                ys.append(kickSize)
#             plt.plot(xs,ys)
#             plt.plot([0,meanDrift+3.01*stdDrift],[thresh,thresh])
#             plt.plot([0,meanDrift+3.01*stdDrift],[thresh*2,thresh*2])
#             plt.show()
    out=out[:-1]+'}'
    out2=out2[:-1]+'}'
    print(out)
    print(out2)

if __name__ == '__main3__':
    #debugConfidenceKick()
    debugAccuracyFreeResponse()
    debugAccuracyInterrogation()
    
if __name__ == '__main2__':
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
    