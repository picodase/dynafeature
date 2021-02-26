# structural.py

'''
structuralPH : A package for structural persistent homology in Python.

Uses algebraic topology tools to create persistence diagrams for protein structures (PDB files).
'''

'''
IMPORT PACKAGES
'''

import Bio.PDB
import numpy as np
import kmapper as km
import scipy
from ripser import ripser
from persim import plot_diagrams
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt

import structurepartitioning as sp
import pdbutils


'''
FUNCTIONS
'''

def runPH(positions:np.array, savefig:bool=False) -> list:
    '''
    Ripser wrapper. Returns the persistence diagram "birth-death" array.
    '''
    try:
        ph = ripser(positions) # compute persistence
        diagrams = ph['dgms'] # extract the persistence diagram
    except:
        return
    
    if savefig:
        plot_diagrams(diagrams, show=False) # plot the diagrams
        plt.savefig("persistent_homology.png") # save the figure
    return diagrams

def visKMapper(data: np.array, id: str):
    '''
    Bundles the functions used to calculate a KMapper visualization. Exports it as an .html object.
    '''

    mapper = km.KeplerMapper(verbose=1)     # init
    projected_data = mapper.fit_transform(data, projection=[0,1]) # fit, transform data to X-Y axis
    graph = mapper.map(projected_data, data,)   # Create dictionary called 'graph' with nodes, edges and meta-information
    mapper.visualize(graph,
                    path_html="make_circles_keplermapper_output"+id+".html",
                    title="make_circles(n_samples=5000, noise=0.03, factor=0.3)")
    return

