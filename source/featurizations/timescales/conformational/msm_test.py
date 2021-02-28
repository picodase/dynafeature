# msm_test.py

from .msmrunner import MSMRunner

# import PDB and trajectories as filepath strings
pdb = 'data/md_0_1_298.pdb'
traj = 'data/md_0_1_298_reduced.xtc'

msmr = MSMRunner(pdb, traj)

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


