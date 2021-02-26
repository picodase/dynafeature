
'''


'''


'''
IMPORTS
'''
import numpy as np
import pandas as pd
from prody import *

'''
FUNCTIONS
'''

def extractModesCoords(filename: str = 'ubi_pca.nmd') -> tuple:
    '''
    Extracts normal modes and atom positions from a nmd file into a tuple of np.ndarrays.
    '''
    # extract atom positions 
    df = pd.read_csv(filename, skiprows=8, sep=" ")
    cols = df.columns.values.tolist()
    for i in list(range(len(cols))):    # for the Ubiquitin test
        if cols[i] == '21.758.1':
            cols[i] = '21.758'
    coords = np.array([float(i) for i in cols[1:]])
    coords = coords.reshape((-1,3))     # reshape
    
    # extract coords and modes
    mode_df = pd.read_csv(filename, skiprows=8, sep=" ").transpose()['mode']
    mode_df.drop(mode_df.index[[0]])
    # make modes numpy array
    modes = mode_df.to_numpy()
    intensities = modes[0]
    modes = modes[1:]
    modes = modes.transpose().reshape((6,-1,3)) # reshape

    return (coords, modes)

def computeNMDfile(filename: str = '2k39'):
    '''
    Runs a NMA analysis on a multi-frame object. Here, we start with an NMR structure.
    '''
    # obtain PCA dataset as a .nmd file
    ubi = parsePDB(filename, subset='calpha')
    ubi_selection = ubi.select('resnum < 71')
    ubi_ensemble = Ensemble(ubi_selection)
    ubi_ensemble.iterpose()
    pca = PCA('ubiquitin')
    pca.buildCovariance(ubi_ensemble)
    pca.calcModes()
    writeNMD('ubi_pca.nmd', pca[:6], ubi_selection)
    return

