

# %%
# imports

import pdbutils
import elemSpecificPH

# %%
# get the structure
struc = pdbutils.createStruct("6VXX.pdb")
# %%
# separate it by its elements
#coords = elemSpecificPH.getPositionArray(struc)

# %%
strucarray = elemSpecificPH.struc2elemPosArrays(struc)
# %%
#a = strucarray[1]

# %%
# for each substructure (element selection) to analyze,
for i in list(range(len(strucarray))):
    # PERSISTENT HOMOLOGY
    elemSpecificPH.printPHDiagram(strucarray[i+1], str(i+1)) # print the persistent homology diagram
    
    # KMAPPER
    #elemSpecificPH.visKMapper(strucarray[i+1], str(i+1))
    # need to resolve some errors...

# %%

