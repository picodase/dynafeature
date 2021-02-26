
# %%
import numpy as np
from ripser import ripser
from persim import plot_diagrams

data = np.random.random((100,2))
diagrams = ripser(data)['dgms']
plot_diagrams(diagrams, show=True)

# %%
cocycles = ripser(data, do_cocycles=True)['cocycles']
cocycles

# %%
ret = ripser(data, do_cocycles=True)
#ret

dgms = ret['dgms']
cocycles = ret['cocycles']
num_edges = ret['num_edges']
dperm2all = ret['dperm2all']
idx_perm = ret['idx_perm']

# %%
dgms
# %%
cocycles

# %%
