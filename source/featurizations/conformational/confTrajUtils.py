
'''
IMPORTS
'''

import Bio.PDB
import numpy as np

import pyEMMA as pe

import structureUtils as su

'''
CONSTANTS
'''



'''
FUNCTIONS
'''

class Conformation():
    def __init__(self, struc:su.Structure):
        self.structure = struc
        self.partitions = ap.AtomPartition(self.structure.biostructure)

class ConformationPair():
    def __init__(self, conf1:Conformation, conf2:Conformation):
        self.c1 = conf1
        self.c2 = conf2
        self.wdistPH = 0
        
class ConformationSet():
    def __init__(self, struclist:list):
        self.conflist = struclist
        

class Trajectory():
    def __init__(self, trajfile:str, topofile:str):
        self.traj = 0
        
def extractConformations(trajfile:str, topofile:str) -> list:
    '''
    DESCRIPTION
        Runs a Markov State Model to extract conformations in an MD trajectory, and returns
        i)      the list of PDB filenames corresponding to conformations,
        ii)     summary plots for the Markov State Model in CWD,

    ARGUMENTS

        trajfile:str    String filename of the trajectory.
        topofile:str    String filename of the topology.
    '''
    confs = []
    
    # Run MSM functions here on traj.
    return confs

