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
    LOGODDS_PERFECT = 8
    LOGODDS_COMBO = 9
    LOGODDS_DISTRIBUTION = 10
    LOGODDS_APPROX = 11
    BELIEF = 12
    
class UpdateModel(Enum):
    NO_UPDATE = 0
    THETA_UPDATE = 1
    THETA_NORM = 2
    OPTIMAL = 3
    BELIEF_UP = 4
    FINITE_TIME = 5

class NetworkType(Enum):
    FULLY_CONNECTED = 0
    ERSOS_RENYI = 1
    BARABASI_ALBERT = 2
    SPACE = 3
    SOFT_RGG = 4
    FROM_FILE = 5
    FROM_FILE_FIXED_COMM = 6
    RGG_FIXED_DEGREE = 7
    
class DriftDistribution(Enum):
    UNIFORM = 0
    NORMAL = 1
    FROM_ACCURACY = 2

class AgentType(Enum):
    SIMPLE = 0
    DDM = 1
    

    