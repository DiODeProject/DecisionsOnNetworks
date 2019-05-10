'''
Created on 27 Jan 2017

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''
from enum import Enum

class DecisionModel(Enum):
    CONFIDENCE = 1
    MAJORITY_RAND = 2
    BELIEF = 12
    
class UpdateModel(Enum):
    NO_UPDATE = 0
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
    
    