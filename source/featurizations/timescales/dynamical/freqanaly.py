

# %% imports

import igraph
import pandas as pd
import matplotlib.pyplot as plt
import cairo

# %% 
## conversion to network

df = pd.read_csv("data/contfreqs_mdsim_07_10_10_13_49_07.tsv", sep="\t", skiprows=2, names=["res1", "res2", "freq"])

# %%
# plot histogram

# %% 
# create network from this data

g = igraph.Graph.TupleList(df.values.tolist(), weights=True)

# %%
g.summary()

# %%
igraph.drawing.plot(g, edge_curved=True)

# %%
thrmot = g.motifs_randesu(size=3)
fourmot = g.motifs_randesu(size=4)
# %%
thrmot