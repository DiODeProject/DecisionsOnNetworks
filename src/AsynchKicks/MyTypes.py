'''
Created on 15 Mar 2019

@author: Andreagiovanni Reina.
University of Sheffield, UK.
'''
from enum import Enum
  
class UpdateModel(Enum):
    NO_UPDATE = 0
    THRESH_KICK = 1
    CONF_KICK = 2

class NetworkType(Enum):
    FULLY_CONNECTED = 0
    ERSOS_RENYI = 1
    BARABASI_ALBERT = 2
    SPACE = 3
    SOFT_RGG = 4
    FROM_FILE = 5
    FROM_FILE_FIXED_COMM = 6
    
class DriftDistribution(Enum):
    UNIFORM = 0
    NORMAL = 1
    FROM_ACCURACY = 2

class AgentType(Enum):
    SIMPLE = 0
    DDM = 1
    

    