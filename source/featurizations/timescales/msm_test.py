# msm_test.py

# %%
import msmrunner

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
#import mdshare
import pyemma
from pyemma.util.contexts import settings
# %%
# import PDB and trajectories as filepath strings
pdb = 'data/md_0_1_298.pdb'
files = 'data/md_0_1_298_reduced.xtc'
# %%
labels, torsions_feat, torsions_data, positions_feat, positions_data, distances_feat, distances_data = msm_runner.featurize()
# %%
# feature selection
dim = 10
dim = 1

msmrunner.figVAMP2(torsions_data, positions_data, distances_data, dim=dim)
# %%
# VAMP2 scoring
lags = [1, 2, 5, 10, 20]
dims = [i + 1 for i in range(10)]
dims = [1]

msmrunner.scoringVAMP2(torsions_data, dim=dims, lag=lags)
# %%
tica = pyemma.coordinates.tica(torsions_data, lag=5)
tica_output = tica.get_output()
tica_concatenated = np.concatenate(tica_output)
# %%
# visualize TICA components
'''
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
pyemma.plots.plot_feature_histograms(
    tica_concatenated,
    ax=axes[0],
    #feature_labels=['IC1', 'IC2', 'IC3', 'IC4'],
    ylog=True)
pyemma.plots.plot_density(*tica_concatenated[:, :2].T, ax=axes[1], logscale=True)
axes[1].set_xlabel('IC 1')
axes[1].set_ylabel('IC 2')
fig.tight_layout()

fig.savefig("TICA_components.png")

'''
# %%
msmrunner.metastability(tica_output)
# %%
# discretization
n_clustercenters = [5, 10, 30, 75]
msmrunner.discretization(tica_output, n_clustercenters)
# %%
# find k!
cluster, dtrajs_concatenated = msmrunner.findk(tica_output, k)
# %%
msmrunner.findTICA2(tica_concatenated, cluster)
# %%
# implied timescales
its = pyemma.msm.its(cluster.dtrajs, lags=50, nits=10, errors='bayes')
pyemma.plots.plot_implied_timescales(its, units='ns', dt=0.1)
# %%
msm = pyemma.msm.bayesian_markov_model(cluster.dtrajs, lag=5, dt_traj='0.1 ns')
print('fraction of states used = {:.2f}'.format(msm.active_state_fraction))
print('fraction of counts used = {:.2f}'.format(msm.active_count_fraction))
# %%
# chapman-kolmogorov
nstates = 5

msmrunner.cktester(msm)
# %%
# spectral analysis
nits = 15

msmrunner.spectralAnalysis(msm)
# %%
# analyze stationary distribution over first PCs
msmrunner.stationaryDistrib(tica_concatenated,msm,dtrajs_concatenated)
# %%
# eigenvectors
msmrunner.principalEigvecs(msm, tica_concatenated, dtrajs_concatenated)
# %%
# PCCA
msm.pcca(nstates)
# %% 
# first 2 TICA dims
msmrunner.first2TICAdims(msm,tica_concatenated,dtrajs_concatenated)
# %%
msmrunner.metastability2(msm, dtrajs_concatenated, tica_concatenated, nstates)
# %%
# write trajs to disk
msmrunner.writeTrajsToDisk(msm, torsions_feat)
# %%
# free energy of the states
msmrunner.plotFreeEnergyDiag(msm)
# %%
# extract MFPT
df = msmrunner.extractMFPT(msm, 'ns')
df.to_csv("MFPT.csv")
# %%
# 
A = msm.metastable_sets[0]
B = np.concatenate(msm.metastable_sets[1:])
print('MFPT 1 -> other: ({:6.1f} ± {:5.1f}) ns'.format(
    msm.sample_mean('mfpt', A, B), msm.sample_std('mfpt', A, B)))
print('MFPT other -> 1: ({:.1f} ± {:5.1f}) ns'.format(
    msm.sample_mean('mfpt', B, A), msm.sample_std('mfpt', B, A)))
# %%
# transition path theory
start, final = 1, 3
A = msm.metastable_sets[start]
B = msm.metastable_sets[final]
flux = pyemma.msm.tpt(msm, A, B)

cg, cgflux = flux.coarse_grain(msm.metastable_sets)
# %%
# Project onto first two TICA components
fig, ax = plt.subplots(figsize=(5, 4))

pyemma.plots.plot_contour(
    *tica_concatenated[:, :2].T,
    flux.committor[dtrajs_concatenated],
    cmap='brg',
    ax=ax,
    mask=True,
    cbar_label=r'committor $\mathcal{S}_%d \to \mathcal{S}_%d$' % (
        start + 1, final + 1))
fig.tight_layout()
fig.savefig("onto_TICA.png")