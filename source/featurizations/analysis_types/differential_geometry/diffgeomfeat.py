

'''
This package includes modification functions to be used in persistent homology computations.

'''

'''
IMPORTS
'''
import numpy as np
import Bio.PDB
import math

'''
FUNCTIONS
'''

### CORRELATION KERNELS

def generalizedExp(r_i: float, r_j: float, η_kk: float, κ: float) -> float:
    '''
    Computes the generalized exponential function for a norm ||r_i - r_j||.
    '''
    return (math.exp(-exp(r_i-r_j)/η_kk)^κ)

def generalizedLorentz(r_i: float, r_j: float, η_kk: float, ν: float) -> float:
    '''
    Computes the generalized Lorentz function for a norm ||r_i - r_j||.
    '''
    return (1 + abs(r_i-r_j)/η_kk)^ν)^-1

def elementInteractiveDensity():
    