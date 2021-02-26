
'''
This package is used to interface with pyemma runner functions to simplify analyses.
'''


'''
IMPORTS
'''

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
#import mdshare
import pyemma
from pyemma.util.contexts import settings

'''
FUNCTIONS
'''

def featurize(pdb: str) -> tuple:
    '''
    Performs featurizations for MSM construction.
    '''

    torsions_feat = pyemma.coordinates.featurizer(pdb)
    torsions_feat.add_backbone_torsions(cossin=True, periodic=False)
    torsions_data = pyemma.coordinates.load(files, features=torsions_feat)
    labels = ['backbone\ntorsions']

    positions_feat = pyemma.coordinates.featurizer(pdb)
    positions_feat.add_selection(positions_feat.select_Backbone())
    positions_data = pyemma.coordinates.load(files, features=positions_feat)
    labels += ['backbone atom\npositions']

    distances_feat = pyemma.coordinates.featurizer(pdb)
    distances_feat.add_distances(
        distances_feat.pairs(distances_feat.select_Backbone(), excluded_neighbors=2), periodic=False)
    distances_data = pyemma.coordinates.load(files, features=distances_feat)
    labels += ['backbone atom\ndistances']

    return (labels, torsions_feat, torsions_data, positions_feat, positions_data, distances_feat, distances_data)

def score_cv(data, dim, lag, number_of_splits=10, validation_fraction=0.5):
    """Compute a cross-validated VAMP2 score.
    
    We randomly split the list of independent trajectories into
    a training and a validation set, compute the VAMP2 score,
    and repeat this process several times.
    
    Parameters
    ----------
    data : list of numpy.ndarrays
        The input data.
    dim : int
        Number of processes to score; equivalent to the dimension
        after projecting the data with VAMP2.
    lag : int
        Lag time for the VAMP2 scoring.
    number_of_splits : int, optional, default=10
        How often do we repeat the splitting and score calculation.
    validation_fraction : int, optional, default=0.5
        Fraction of trajectories which should go into the validation
        set during a split.
    """
    # we temporarily suppress very short-lived progress bars
    with pyemma.util.contexts.settings(show_progress_bars=False):
        nval = int(len(data) * validation_fraction)
        scores = np.zeros(number_of_splits)
        for n in range(number_of_splits):
            ival = np.random.choice(len(data), size=nval, replace=False)
            vamp = pyemma.coordinates.vamp(
                [d for i, d in enumerate(data) if i not in ival], lag=lag, dim=dim)
            scores[n] = vamp.score([d for i, d in enumerate(data) if i in ival])
    return scores

def figVAMP2(torsions_data, positions_data, distances_data, dim):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3), sharey=True)
    for ax, lag in zip(axes.flat, [5, 10, 20]):
        torsions_scores = score_cv(torsions_data, lag=lag, dim=dim)
        scores = [torsions_scores.mean()]
        errors = [torsions_scores.std()]
        positions_scores = score_cv(positions_data, lag=lag, dim=dim)
        scores += [positions_scores.mean()]
        errors += [positions_scores.std()]
        distances_scores = score_cv(distances_data, lag=lag, dim=dim)
        scores += [distances_scores.mean()]
        errors += [distances_scores.std()]
        ax.bar(labels, scores, yerr=errors, color=['C0', 'C1', 'C2'])
        ax.set_title(r'lag time $\tau$={:.1f}ns'.format(lag * 0.1))
        if lag == 5:
            # save for later
            vamp_bars_plot = dict(
                labels=labels, scores=scores, errors=errors, dim=dim, lag=lag)
    axes[0].set_ylabel('VAMP2 score')
    fig.tight_layout()
    fig.savefig("VAMP2.png")
    return

def scoringVAMP2(torsions_data, dims, lags):
    fig, ax = plt.subplots()
    for i, lag in enumerate(lags):
        scores_ = np.array([score_cv(torsions_data, dim, lag)
                            for dim in dims])
        scores = np.mean(scores_, axis=1)
        errors = np.std(scores_, axis=1, ddof=1)
        color = 'C{}'.format(i)
        ax.fill_between(dims, scores - errors, scores + errors, alpha=0.3, facecolor=color)
        ax.plot(dims, scores, '--o', color=color, label='lag={:.1f}ns'.format(lag * 0.1))
    ax.legend()
    ax.set_xlabel('number of dimensions')
    ax.set_ylabel('VAMP2 score')
    fig.tight_layout()

    fig.savefig("VAMP2_2.png")
    return

def metastability(tica_output):
    fig, axes = plt.subplots(4, 1, figsize=(12, 5), sharex=True)
    x = 0.1 * np.arange(tica_output[0].shape[0])
    for i, (ax, tic) in enumerate(zip(axes.flat, tica_output[0].T)):
        ax.plot(x, tic)
        ax.set_ylabel('IC {}'.format(i + 1))
    axes[-1].set_xlabel('time / ns')
    fig.tight_layout()

    fig.savefig("metastability.png")
    return

def discretization(tica_output, n_clustercenters):
    scores = np.zeros((len(n_clustercenters), 5))
    for n, k in enumerate(n_clustercenters):
        for m in range(5):
            with pyemma.util.contexts.settings(show_progress_bars=False):
                _cl = pyemma.coordinates.cluster_kmeans(
                    tica_output, k=k, max_iter=50, stride=50)
                _msm = pyemma.msm.estimate_markov_model(_cl.dtrajs, 5)
                scores[n, m] = _msm.score_cv(
                    _cl.dtrajs, n=1, score_method='VAMP2', score_k=min(10, k))

    fig, ax = plt.subplots()
    lower, upper = pyemma.util.statistics.confidence_interval(scores.T.tolist(), conf=0.9)
    ax.fill_between(n_clustercenters, lower, upper, alpha=0.3)
    ax.plot(n_clustercenters, np.mean(scores, axis=1), '-o')
    ax.semilogx()
    ax.set_xlabel('number of cluster centers')
    ax.set_ylabel('VAMP-2 score')
    fig.tight_layout()

    fig.savefig("discretization.png")
    return

def findk(tica_output, k): 
    cluster = pyemma.coordinates.cluster_kmeans(
    tica_output, k=k, max_iter=50, stride=10)
    dtrajs_concatenated = np.concatenate(cluster.dtrajs)
    return (cluster, dtrajs_concatenated)

def findTICA2(tica_concatenated, cluster):
    ig, ax = plt.subplots(figsize=(4, 4))
    pyemma.plots.plot_density(
        *tica_concatenated[:, :2].T, ax=ax, cbar=False, alpha=0.3)
    ax.scatter(*cluster.clustercenters[:, :2].T, s=5, c='C1')
    ax.set_xlabel('IC 1')
    ax.set_ylabel('IC 2')
    fig.tight_layout()

    fig.savefig("TICA_2.png")
    return 

def spectralAnalysis(msm):
    timescales_mean = msm.sample_mean('timescales', k=nits)
    timescales_std = msm.sample_std('timescales', k=nits)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    
    def its_separation_err(ts, ts_err):
        """
        Error propagation from ITS standard deviation to timescale separation.
        """
        return ts[:-1] / ts[1:] * np.sqrt(
            (ts_err[:-1] / ts[:-1])**2 + (ts_err[1:] / ts[1:])**2)

    axes[0].errorbar(
        range(1, nits + 1),
        timescales_mean, 
        yerr=timescales_std, 
        fmt='.', markersize=10)
    axes[1].errorbar(
        range(1, nits),
        timescales_mean[:-1] / timescales_mean[1:], 
        yerr=its_separation_err(
            timescales_mean, 
            timescales_std), 
        fmt='.', 
        markersize=10,
        color='C0')

    for i, ax in enumerate(axes):
        ax.set_xticks(range(1, nits + 1))
        ax.grid(True, axis='x', linestyle=':')
        
    axes[0].axhline(msm.lag * 0.1, lw=1.5, color='k')
    axes[0].axhspan(0, msm.lag * 0.1, alpha=0.3, color='k')
    axes[0].set_xlabel('implied timescale index')
    axes[0].set_ylabel('implied timescales / ns')
    axes[1].set_xticks(range(1, nits))
    axes[1].set_xticklabels(
        ["{:d}/{:d}".format(k, k + 1) for k in range(1, nits)],
        rotation=45)
    axes[1].set_xlabel('implied timescale indices')
    axes[1].set_ylabel('timescale separation')
    fig.tight_layout()

    fig.savefig("spectral_analysis.png")
    return

def cktester(msm):
    cktest = msm.cktest(nstates, mlags=6)
    pyemma.plots.plot_cktest(cktest, dt=0.1, units='ns');
    fig.savefig("cktest.png")
    return

def stationaryDistrib(tica_concatenated, msm, dtrajs_concatenated):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
    pyemma.plots.plot_contour(
        *tica_concatenated[:, :2].T,
        msm.pi[dtrajs_concatenated],
        ax=axes[0],
        mask=True,
        cbar_label='stationary distribution')
    pyemma.plots.plot_free_energy(
        *tica_concatenated[:, :2].T,
        weights=np.concatenate(msm.trajectory_weights()),
        ax=axes[1],
        legacy=False)
    for ax in axes.flat:
        ax.set_xlabel('IC 1')
    axes[0].set_ylabel('IC 2')
    axes[0].set_title('Stationary distribution', fontweight='bold')
    axes[1].set_title('Reweighted free energy surface', fontweight='bold')
    fig.tight_layout()

    fig.savefig("stationary_distrib.png")
    return

def principalEigvecs(msm, tica_concatenated, dtrajs_concatenated):
    eigvec = msm.eigenvectors_right()
    print('The first eigenvector is one: {} (min={}, max={})'.format(
        np.allclose(eigvec[:, 0], 1, atol=1e-15), eigvec[:, 0].min(), eigvec[:, 0].max()))

    fig, axes = plt.subplots(1, 4, figsize=(15, 3), sharex=True, sharey=True)
    for i, ax in enumerate(axes.flat):
        pyemma.plots.plot_contour(
            *tica_concatenated[:, :2].T,
            eigvec[dtrajs_concatenated, i + 1],
            ax=ax,
            cmap='PiYG',
            cbar_label='{}. right eigenvector'.format(i + 2),
            mask=True)
        ax.set_xlabel('IC 1')
    axes[0].set_ylabel('IC 2')
    fig.tight_layout()

    fig.savefig("principal_eigvecs.png")
    return

def first2TICAdims(msm, tica_concatenated,dtrajs_concatenated):
    fig, axes = plt.subplots(1, 5, figsize=(15, 3), sharex=True, sharey=True)
    for i, ax in enumerate(axes.flat):
        pyemma.plots.plot_contour(
            *tica_concatenated[:, :2].T,
            msm.metastable_distributions[i][dtrajs_concatenated],
            ax=ax,
            cmap='afmhot_r', 
            mask=True,
            cbar_label='metastable distribution {}'.format(i + 1))
        ax.set_xlabel('IC 1')
    axes[0].set_ylabel('IC 2')
    fig.tight_layout()

    fig.savefig("first_2_TICA_dims.png")
    return

def metastability2(msm, dtrajs_concatenated, tica_concatenated, nstates):
    metastable_traj = msm.metastable_assignments[dtrajs_concatenated]

    fig, ax = plt.subplots(figsize=(5, 4))
    _, _, misc = pyemma.plots.plot_state_map(
        *tica_concatenated[:, :2].T, metastable_traj, ax=ax)
    ax.set_xlabel('IC 1')
    ax.set_ylabel('IC 2')
    misc['cbar'].set_ticklabels([r'$\mathcal{S}_%d$' % (i + 1)
                                for i in range(nstates)])
    fig.tight_layout()

    fig.savefig("metastability2.png")

def writeTrajsToDisk(msm, torsions_feat):
    pcca_samples = msm.sample_by_distributions(msm.metastable_distributions, 10)
    torsions_source = pyemma.coordinates.source(files, features=torsions_feat)
    pyemma.coordinates.save_trajs(
        torsions_source,
        pcca_samples,
        outfiles=['pcca{}_10samples.pdb'.format(n + 1)
                for n in range(msm.n_metastable)])
    return

def plotFreeEnergyDiag(msm):
    print('state\tπ\t\tG/kT')
    for i, s in enumerate(msm.metastable_sets):
        p = msm.pi[s].sum()
        print('{}\t{:f}\t{:f}'.format(i + 1, p, -np.log(p)))
    return

def extractMFPT(msm, units: str) -> pd.DataFrame:
    from itertools import product

    mfpt = np.zeros((nstates, nstates))
    for i, j in product(range(nstates), repeat=2):
        mfpt[i, j] = msm.mfpt(
            msm.metastable_sets[i],
            msm.metastable_sets[j])

    from pandas import DataFrame
    print('MFPT / '+units+':')
    df = DataFrame(np.round(mfpt, decimals=2), index=range(1, nstates + 1), columns=range(1, nstates + 1))
    
    return df