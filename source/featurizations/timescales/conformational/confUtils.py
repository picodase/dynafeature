
'''
IMPORTS
'''

import numpy as np

from ..structural.structureUtils import Structure

'''
CONSTANTS
'''



'''
FUNCTIONS
'''

class Conformation(Structure):
    def __init__(self):
        self.partitions = ap.AtomPartition(self.structure.biostructure)

class ConformationPair():
    def __init__(self, conf1:Conformation, conf2:Conformation):
        self.c1 = conf1
        self.c2 = conf2

class ConformationSet():
    def __init__(self, struclist:list):
        self.conflist = struclist

