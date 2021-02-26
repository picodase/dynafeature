
import pyemma
import Bio.PDB

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

