# %% 
# imports
import nmafeat
from prody import *
import numpy as np
# %%
nmafeat.computeNMDfile()
(coords, modes) = nmafeat.extractModesCoords()
print(f"Values range from {coords.min()} to {coords.max()}")
# %%
dmin = coords.min() - 2
dmax = coords.max() + 2
a = np.mgrid[dmin:dmax:101j,dmin:dmax:101j,dmin:dmax:101j]
a = a.reshape(3,1030301).transpose()
# %%
import scipy.interpolate.ndgriddata as ndgriddata
mode1_interp = ndgriddata.griddata(coords, modes[0], a)     # interpolate the data onto the grid
# %%
# Eureka! Maybe it worked!: 
# is not nan
mode1_interp[~np.isnan(mode1_interp)].size
# %%
# is nan
mode1_interp[np.isnan(mode1_interp)].size
# %%
δ = np.gradient(mode1_interp)
# %%
δ[0][~np.isnan(δ[0])].size
# %%
m1i_orig = mode1_interp.reshape(3,101,101,101)
# %%
m1i_grad = np.gradient(m1i_orig)
# %%
m1i_grad[0].shape
# %%
# find where the gradient is zero
np.where(m1i_grad[0] == 0)
# %%
