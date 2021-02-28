'''
IMPORTS
'''
import Bio.PDB
import numpy as np
#import pdbutils

from .structureUtils import Structure

'''
CONSTANTS
'''

carbontypes = {"C","CA","CB","CG","CE"}
hydrophobics = {"PHE","ALA","LEU","MET","ILE","TRP","PRO","VAL"}
polarunchargeds = {"CYS","GLY","GLN","ASN","SER","TYR","THR"}
acidics = {"ASP","GLU"}
basics = {"HIS","ARG","LYS"}
larges = {"ARG","LYS","TRP","MET","GLU","GLN","TYR","PHE"}
smalls = {"HIS","ILE","LEU","ALA","PRO","VAL","CYS","GLY","SER","ASN","THR","ASP"}
specials = {"GLY", "CYS", "PRO"}

'''
FUNCTIONS
'''

def getAllAtomPositions(struc:Bio.PDB.Structure.Structure) -> np.ndarray:
    '''
    Obtains a position array for the 3D-coordinates of all atoms in the provided structure.
    '''
    coords = []   # construct a (3 x natoms) array
    # take the structure object
    # faster to do it this way, see:
    # https://stackoverflow.com/questions/3881453/numpy-add-row-to-array

    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                for atom in res:    # for each atom,
                    coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray

class StructurePartition():
    def __init__(self, struc:Structure):
        self.structure = struc

class MacromoleculePartition(StructurePartition):
    def __init__(self, struc:Structure, it:str="only_protein"):
        self.interactype = it   # Interaction type
        self.structure = struc
        
class onlyLigand(MacromoleculePartition):
    def __init__(self, struc:Structure):

class onlyProtein(MacromoleculePartition):
    def __init__(self, struc):
        self.bindingpocket = 0
        self.distant = 0

class interactionsLigandProtein(MacromoleculePartition):
    def __init__(self, struc):
        self.firstshell = 0
        self.secondshell = 0
        self.distant = 0
            
class AtomPartition(StructurePartition):
    def sepEN(self) -> np.ndarray:
        '''
        Obtains a position array for the 3D-coordinates of each atom in the provided structure.
        '''
        coords = []   # construct a (3 x natoms) array
        for model in self.structure: # for each model,
            for chain in model:    # for each chain,
                for res in chain:   # for each residue,
                    for atom in res:    # for each atom,
                        if atom.get_id() == "N": # if the atom is an carbon,
                            coords.append(atom.get_coord()) # get the coordinates, append to a new array
        return np.array(coords)    # convert to np.ndarray

    def sepEC(self) -> np.ndarray:
        '''
        Obtains a position array for the 3D-coordinates of each atom in the provided structure.
        '''
        coords = []   # construct a (3 x natoms) array
        for model in self.structure: # for each model,
            for chain in model:    # for each chain,
                for res in chain:   # for each residue,
                    for atom in res:    # for each atom,
                        if atom.get_id() in carbontypes: # if the atom is a carbon,
                            coords.append(atom.get_coord()) # get the coordinates, append to a new array
        return np.array(coords)    # convert to np.ndarray

    def sepEO(self) -> np.ndarray:
        '''
        Obtains a position array for the 3D-coordinates of each atom in the provided structure.
        '''
        coords = []   # construct a (3 x natoms) array
        for model in self.structure: # for each model,
            for chain in model:    # for each chain,
                for res in chain:   # for each residue,
                    for atom in res:    # for each atom,
                        if atom.get_id() == "O": # if the atom is an oxygen,
                            coords.append(atom.get_coord()) # get the coordinates, append to a new array
        return np.array(coords)    # convert to np.ndarray

    def sepES(self) -> np.ndarray:
        '''
        Obtains a position array for the 3D-coordinates of each atom in the provided structure.
        '''
        coords = []   # construct a (3 x natoms) array
        for model in self.structure: # for each model,
            for chain in model:    # for each chain,
                for res in chain:   # for each residue,
                    for atom in res:    # for each atom,
                        if atom.get_id() == "S": # if the atom is a sulfur,
                            coords.append(atom.get_coord()) # get the coordinates, append to a new array
        return np.array(coords)    # convert to np.ndarray

    ###### AMINO ACID TYPES ######

    def sepRHψ(self) -> np.ndarray:
        '''
        Obtains a position array for the 3D-coordinates of each atom in the provided structure.
        '''
        coords = []   # construct a (3 x natoms) array
        for model in self.structure: # for each model,
            for chain in model:    # for each chain,
                for res in chain:   # for each residue,
                    if res.get_resname() in hydrophobics: # if residue is in the class of hydrophobic aa's,
                        for atom in res:    # for each atom,
                            coords.append(atom.get_coord()) # get the coordinates, append to a new array
        return np.array(coords)    # convert to np.ndarray

    def sepRδ0(self) -> np.ndarray:
        '''
        Obtains a position array for the 3D-coordinates of each atom in the provided structure.
        '''
        coords = []   # construct a (3 x natoms) array
        for model in self.structure: # for each model,
            for chain in model:    # for each chain,
                for res in chain:   # for each residue,
                    if res.get_resname() in polarunchargeds:  # if residue is in the class of polar aa's,
                        for atom in res:    # for each atom,
                            coords.append(atom.get_coord()) # get the coordinates, append to a new array
        return np.array(coords)    # convert to np.ndarray

    def sepRδa(self) -> np.ndarray:
        '''
        Obtains a position array for the 3D-coordinates of each atom in the provided structure.
        '''
        coords = []   # construct a (3 x natoms) array
        for model in self.structure: # for each model,
            for chain in model:    # for each chain,
                for res in chain:   # for each residue,
                    if res.get_resname() in acidics:  # if residue is in the class of polar aa's,
                        for atom in res:    # for each atom,
                            coords.append(atom.get_coord()) # get the coordinates, append to a new array
        return np.array(coords)    # convert to np.ndarray

    def sepRδb(self) -> np.ndarray:
        '''
        Obtains a position array for the 3D-coordinates of each atom in the provided structure.
        '''
        coords = []   # construct a (3 x natoms) array
        for model in self.structure: # for each model,
            for chain in model:    # for each chain,
                for res in chain:   # for each residue,
                    if res.get_resname() in basics:   # if residue is in the class of polar aa's,
                        for atom in res:    # for each atom,
                            coords.append(atom.get_coord()) # get the coordinates, append to a new array
        return np.array(coords)    # convert to np.ndarray

    def sepRl(self) -> np.ndarray:
        '''
        Obtains a position array for the 3D-coordinates of each atom in the provided structure.
        '''
        coords = []   # construct a (3 x natoms) array
        for model in self.structure: # for each model,
            for chain in model:    # for each chain,
                for res in chain:   # for each residue,
                    if res.get_resname() in larges:   # if residue is in the class of large aa's,
                        for atom in res:    # for each atom,
                            coords.append(atom.get_coord()) # get the coordinates, append to a new array
        return np.array(coords)    # convert to np.ndarray

    def sepRs(self) -> np.ndarray:
        '''
        Obtains a position array for the 3D-coordinates of each atom in the provided structure.
        '''
        coords = []   # construct a (3 x natoms) array
        for model in self.structure: # for each model,
            for chain in model:    # for each chain,
                for res in chain:   # for each residue,
                    if res.get_resname() in smalls: # if residue is in the class of small aa's,
                        for atom in res:    # for each atom,
                            coords.append(atom.get_coord()) # get the coordinates, append to a new array
        return np.array(coords)    # convert to np.ndarray


    ###### TYPE ######

    def sepTs(self) -> np.ndarray:
        '''
        Obtains a position array for the 3D-coordinates of each atom in the provided structure.
        '''
        coords = []   # construct a (3 x natoms) array
        for model in self.structure: # for each model,
            for chain in model:    # for each chain,
                for res in chain:   # for each residue,
                    if res.get_resname() in specials: # if residue is in the class of small aa's,
                        for atom in res:    # for each atom,
                            coords.append(atom.get_coord()) # get the coordinates, append to a new array
        return np.array(coords)    # convert to np.ndarray
    
    def __init__(self, s:Structure):

        self.byelemtype = [self.sepEC(), self.sepEN(), self.sepEO(), self.sepES()]
        self.byaatype = [self.sepRHψ(),self.sepRl(),self.sepRs(),self.sepRδ0(),self.sepRδa(),self.sepRδb()]
        self.bysstype = [self.sepTs()]
    
    def getAllSeps(self):
        return [self.byelemtype,self.byaatype,self.bysstype]

###### SECONDARY STRUCTURE ######
'''
def sssphβ(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    #Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res # if residue is in a beta-sheet,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray

# modify the sep function

def sssphα(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    #Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res # if residue is in an alpha-helix,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray


def sssphLoop(struc: Bio.PDB.Structure.Structure) -> np.ndarray:

    #Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res # if residue is in a loop,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray
    '''
