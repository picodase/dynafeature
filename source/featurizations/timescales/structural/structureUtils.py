

'''
IMPORTS
'''

import Bio.PDB

from Bio.PDB.PDBParser import PDBParser

import numpy as np

'''
FUNCTIONS
'''

def Structure():
    def __init__(self, filename:str):
        parser = PDBParser(PERMISSIVE=1)
        self.biostruct = parser.get_structure(filename.replace(".pdb", ""), filename)
        
