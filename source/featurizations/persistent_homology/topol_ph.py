
'''
IMPORTS
'''
import pdbutils
import elemSpecificPH
import scipy
from ripser import ripser
from persim import plot_diagrams
import matplotlib.pyplot as plt

'''
FUNCTIONS
'''


'''
MAIN
'''

### SET UP THE STRUCTURE OBJECT -- A NUMPY ARRAY WITH ATOM POSITIONS
struc = pdbutils.createStruct("6VXX.pdb")               # get the structure
coords = elemSpecificPH.getPositionArray(struc)         # separate it by its elements
strucarray = elemSpecificPH.struc2elemPosArrays(struc)  # construct the structure array
a = strucarray[1]                                       # take the first element from the structure array
#a = scipy.spatial.distance.pdist(strucarray)



### RUN PERSISTENT HOMOLOGY COMPUTATIONS
diagrams = ripser(a)['dgms']                            # compute the persistence diagram of the 
plot_diagrams(diagrams, show=False)                     # plot the diagrams
plt.savefig("persistent_homology.png")                  # save the figure
