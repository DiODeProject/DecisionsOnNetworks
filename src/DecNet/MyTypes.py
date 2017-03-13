'''
Created on 27 Jan 2017

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''
from enum import Enum

class DecisionModel(Enum):
    CONFIDENCE = 1
    MAJORITY_RAND = 2
    MAJORITY_BIAS = 3
    MAJORITY_INHIB= 4
    CDCI = 5
    BEST_ACC = 6
    BEST_CONF = 7
    LOGODDS = 8

class NetworkType(Enum):
    FULLY_CONNECTED = 0
    ERSOS_RENYI = 1
    BARABASI_ALBERT = 2
    
class DriftDistribution(Enum):
    UNIFORM = 0
    NORMAL = 1
    FROM_ACCURACY = 2

class AgentType(Enum):
    SIMPLE = 0
    DDM = 1

    