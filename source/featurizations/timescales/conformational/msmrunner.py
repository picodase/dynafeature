
'''
This package is used to interface with pyemma runner functions to simplify analyses.
'''


'''
IMPORTS
'''

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd

import pyemma
from pyemma.util.contexts import settings

'''
FUNCTIONS
'''

class MSMRunner():
    '''
    A class for running Markov State Model analyses--simplified from tutorial code (PyEMMA).
    '''
    def featurize(self):
        '''
        Performs featurizations for MSM construction.
        '''
        torsions_feat = pyemma.coordinates.featurizer(self.topol)
        torsions_feat.add_backbone_torsions(cossin=True, periodic=False)
        torsions_data = pyemma.coordinates.load(self.files, features=torsions_feat)
        labels = ['backbone\ntorsions']

        positions_feat = pyemma.coordinates.featurizer(self.topol)
        positions_feat.add_selection(positions_feat.select_Backbone())
        positions_data = pyemma.coordinates.load(self.files, features=positions_feat)
        labels += ['backbone atom\npositions']

        distances_feat = pyemma.coordinates.featurizer(self.topol)
        distances_feat.add_distances(
            distances_feat.pairs(distances_feat.select_Backbone(), excluded_neighbors=2), periodic=False)
        distances_data = pyemma.coordinates.load(self.files, features=distances_feat)
        labels += ['backbone atom\ndistances']
        self.features = (labels, torsions_feat, torsions_data, positions_feat, positions_data, distances_feat, distances_data)
        return

    def score_cv(self, data, lag, dim:int=10, number_of_splits=10, validation_fraction=0.5):
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
        self.scores = scores
        return

    def figVAMP2(self):
        torsions_data = self.features[2]
        positions_data = self.features[4]
        distances_data = self.features[6]
        fig, axes = plt.subplots(1, 3, figsize=(12, 3), sharey=True)
        for ax, lag in zip(axes.flat, [5, 10, 20]):
            torsions_scores = self.score_cv(torsions_data, lag=lag, dim=self.VAMP2dims)
            scores = [torsions_scores.mean()]
            errors = [torsions_scores.std()]
            positions_scores = self.score_cv(positions_data, lag=lag, dim=self.VAMP2dims)
            scores += [positions_scores.mean()]
            errors += [positions_scores.std()]
            distances_scores = self.score_cv(distances_data, lag=lag, dim=self.VAMP2dims)
            scores += [distances_scores.mean()]
            errors += [distances_scores.std()]
            ax.bar(labels, scores, yerr=errors, color=['C0', 'C1', 'C2'])
            ax.set_title(r'lag time $\tau$={:.1f}ns'.format(lag * 0.1))
            if lag == 5:
                # save for later
                vamp_bars_plot = dict(
                    labels=labels, scores=scores, errors=errors, dim=self.VAMP2dims, lag=lag)
        axes[0].set_ylabel('VAMP2 score')
        fig.tight_layout()
        fig.savefig("VAMP2.png")
        return

    def scoringVAMP2(self):
        fig, ax = plt.subplots()
        for i, lag in enumerate(self.VAMP2lags):
            scores_ = np.array([self.score_cv(self.features[2], dim, lag)
                                for dim in self.dims])
            scores = np.mean(scores_, axis=1)
            errors = np.std(scores_, axis=1, ddof=1)
            color = 'C{}'.format(i)
            ax.fill_between(self.dims, scores - errors, scores + errors, alpha=0.3, facecolor=color)
            ax.plot(self.dims, scores, '--o', color=color, label='lag={:.1f}ns'.format(lag * 0.1))
        ax.legend()
        ax.set_xlabel('number of dimensions')
        ax.set_ylabel('VAMP2 score')
        fig.tight_layout()
        fig.savefig("VAMP2_2.png")
        return

    def metastability(self):
        fig, axes = plt.subplots(4, 1, figsize=(12, 5), sharex=True)
        x = 0.1 * np.arange(self.tica_output[0].shape[0])
        for i, (ax, tic) in enumerate(zip(axes.flat, self.tica_output[0].T)):
            ax.plot(x, tic)
            ax.set_ylabel('IC {}'.format(i + 1))
        axes[-1].set_xlabel('time / ns')
        fig.tight_layout()
        fig.savefig("metastability.png")
        return

    def discretization(self):
        scores = np.zeros((len(self.n_clustercenters), 5))
        for n, k in enumerate(self.n_clustercenters):
            for m in range(5):
                with pyemma.util.contexts.settings(show_progress_bars=False):
                    _cl = pyemma.coordinates.cluster_kmeans(
                        self.tica_output, k=k, max_iter=50, stride=50)
                    _msm = pyemma.msm.estimate_markov_model(_cl.dtrajs, 5)
                    scores[n, m] = _msm.score_cv(
                        _cl.dtrajs, n=1, score_method='VAMP2', score_k=min(10, k))

        fig, ax = plt.subplots()
        lower, upper = pyemma.util.statistics.confidence_interval(scores.T.tolist(), conf=0.9)
        ax.fill_between(self.n_clustercenters, lower, upper, alpha=0.3)
        ax.plot(self.n_clustercenters, np.mean(scores, axis=1), '-o')
        ax.semilogx()
        ax.set_xlabel('number of cluster centers')
        ax.set_ylabel('VAMP-2 score')
        fig.tight_layout()
        fig.savefig("discretization.png")
        return

    def findk(self):
        self.cluster = pyemma.coordinates.cluster_kmeans(
        self.tica_output, k=self.k, max_iter=50, stride=10)
        self.dtrajs_concatenated = np.concatenate(self.cluster.dtrajs)
        return

    def findTICA2(self):
        ig, ax = plt.subplots(figsize=(4, 4))
        pyemma.plots.plot_density(
            *self.tica_concatenated[:, :2].T, ax=ax, cbar=False, alpha=0.3)
        ax.scatter(*self.cluster.clustercenters[:, :2].T, s=5, c='C1')
        ax.set_xlabel('IC 1')
        ax.set_ylabel('IC 2')
        fig.tight_layout()
        fig.savefig("TICA_2.png")
        return

    def spectralAnalysis(self):
        timescales_mean = self.msm.sample_mean('timescales', k=nits)
        timescales_std = self.msm.sample_std('timescales', k=nits)

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

    def cktester(self):
        cktest = self.msm.cktest(self.nstates, mlags=6)
        pyemma.plots.plot_cktest(cktest, dt=0.1, units='ns');
        fig.savefig("cktest.png")
        return

    def stationaryDistrib(self):
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
        pyemma.plots.plot_contour(
            *self.tica_concatenated[:, :2].T,
            self.msm.pi[self.dtrajs_concatenated],
            ax=axes[0],
            mask=True,
            cbar_label='stationary distribution')
        pyemma.plots.plot_free_energy(
            *self.tica_concatenated[:, :2].T,
            weights=np.concatenate(self.msm.trajectory_weights()),
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

    def principalEigvecs(self):
        eigvec = self.msm.eigenvectors_right()
        print('The first eigenvector is one: {} (min={}, max={})'.format(
            np.allclose(eigvec[:, 0], 1, atol=1e-15), eigvec[:, 0].min(), eigvec[:, 0].max()))

        fig, axes = plt.subplots(1, 4, figsize=(15, 3), sharex=True, sharey=True)
        for i, ax in enumerate(axes.flat):
            pyemma.plots.plot_contour(
                *self.tica_concatenated[:, :2].T,
                eigvec[self.dtrajs_concatenated, i + 1],
                ax=ax,
                cmap='PiYG',
                cbar_label='{}. right eigenvector'.format(i + 2),
                mask=True)
            ax.set_xlabel('IC 1')
        axes[0].set_ylabel('IC 2')
        fig.tight_layout()
        fig.savefig("principal_eigvecs.png")
        return

    def first2TICAdims(self):
        fig, axes = plt.subplots(1, 5, figsize=(15, 3), sharex=True, sharey=True)
        for i, ax in enumerate(axes.flat):
            pyemma.plots.plot_contour(
                *self.tica_concatenated[:, :2].T,
                self.msm.metastable_distributions[i][self.dtrajs_concatenated],
                ax=ax,
                cmap='afmhot_r',
                mask=True,
                cbar_label='metastable distribution {}'.format(i + 1))
            ax.set_xlabel('IC 1')
        axes[0].set_ylabel('IC 2')
        fig.tight_layout()
        fig.savefig("first_2_TICA_dims.png")
        return

    def metastability2(self):
        metastable_traj = self.msm.metastable_assignments[self.dtrajs_concatenated]

        fig, ax = plt.subplots(figsize=(5, 4))
        _, _, misc = pyemma.plots.plot_state_map(
            *self.tica_concatenated[:, :2].T, metastable_traj, ax=ax)
        ax.set_xlabel('IC 1')
        ax.set_ylabel('IC 2')
        misc['cbar'].set_ticklabels([r'$\mathcal{S}_%d$' % (i + 1)
                                    for i in range(self.nstates)])
        fig.tight_layout()
        fig.savefig("metastability2.png")
        return

    def writeTrajsToDisk(self):
        pcca_samples = self.msm.sample_by_distributions(self.msm.metastable_distributions, 10)
        torsions_source = pyemma.coordinates.source(self.files, features=self.torsions_feat)
        pyemma.coordinates.save_trajs(
            torsions_source,
            pcca_samples,
            outfiles=['pcca{}_10samples.pdb'.format(n + 1)
                    for n in range(self.msm.n_metastable)])
        return

    def plotFreeEnergyDiag(self):
        print('state\tπ\t\tG/kT')
        for i, s in enumerate(self.msm.metastable_sets):
            p = self.msm.pi[s].sum()
            print('{}\t{:f}\t{:f}'.format(i + 1, p, -np.log(p)))
        return

    def extractMFPT(self, units:str="ns") -> pd.DataFrame:
        from itertools import product

        mfpt = np.zeros((nstates, nstates))
        for i, j in product(range(nstates), repeat=2):
            mfpt[i, j] = self.msm.mfpt(
                self.msm.metastable_sets[i],
                self.msm.metastable_sets[j])

        from pandas import DataFrame
        print('MFPT / '+units+':')
        self.df = DataFrame(np.round(mfpt, decimals=2), index=range(1, nstates + 1), columns=range(1, nstates + 1))
        return

    def __init__(self, top:str, trj:str, dims:int=10):
        self.topol = top
        self.traj = trj
        self.files = [trj]
        self.features = 0
        self.featurize()
        self.dims = 10
        self.dim = 1
        self.VAMP2lags = [1, 2, 5, 10, 20]
        self.VAMP2dims = [i + 1 for i in range(10)]
        self.VAMP2dims = [1]
        self.figVAMP2()
        self.scoringVAMP2(self.features[2], dim=self.VAMP2dims, lag=self.VAMP2lags)
        self.tica = pyemma.coordinates.tica(self.features[2], lag=5)
        self.tica_output = tica.get_output()
        self.tica_concatenated = np.concatenate(self.tica_output)
        self.metastability(self.tica_output)
        self.n_clustercenters = [5, 10, 30, 75]
        self.discretization(self.tica_output, self.n_clustercenters)
        self.k = 10
        self.cluster, self.dtrajs_concatenated = self.findk(self.tica_output, k)
        self.findTICA2(self.tica_concatenated, self.cluster)
        self.its = pyemma.msm.its(cluster.dtrajs, lags=50, nits=10, errors='bayes')
        pyemma.plots.plot_implied_timescales(its, units='ns', dt=0.1)
        self.msm = pyemma.msm.bayesian_markov_model(self.cluster.dtrajs, lag=5, dt_traj='0.1 ns')
        print('fraction of states used = {:.2f}'.format(msm.active_state_fraction))
        print('fraction of counts used = {:.2f}'.format(msm.active_count_fraction))
        self.nstates = 5
        self.cktester(self.msm)
        self.nits = 15
        self.spectralAnalysis(self.msm)
        self.stationaryDistrib(self.tica_concatenated, self.msm, self.dtrajs_concatenated)
        self.principalEigvecs(self.msm, self.tica_concatenated, self.dtrajs_concatenated)
        self.msm.pcca(self.nstates)
        self.first2TICAdims(self.msm, self.tica_concatenated,   self.dtrajs_concatenated)
        self.metastability2(self.msm, self.dtrajs_concatenated, self.tica_concatenated, self.nstates)
        self.writeTrajsToDisk(self.msm, self.features[2])
        self.plotFreeEnergyDiag(msm)
        self.df = self.extractMFPT(msm, 'ns')
        self.df.to_csv("MFPT.csv")
        self.A = msm.metastable_sets[0]
        self.B = np.concatenate(msm.metastable_sets[1:])
        print('MFPT 1 -> other: ({:6.1f} ± {:5.1f}) ns'.format(self/msm.sample_mean('mfpt', self.A, self.B), self.msm.sample_std('mfpt', A, B)))
        print('MFPT other -> 1: ({:.1f} ± {:5.1f}) ns'.format(
        self.msm.sample_mean('mfpt', B, A), self.msm.sample_std('mfpt', B, A)))
        self.start, self.final = 1, 3
        self.A = self.msm.metastable_sets[self.start]
        self.B = self.msm.metastable_sets[self.final]
        self.flux = pyemma.msm.tpt(self.msm, self.A, self.B)
        self.cg, self.cgflux = self.flux.coarse_grain(self.msm.metastable_sets)
        fig, ax = plt.subplots(figsize=(5, 4))
        pyemma.plots.plot_contour(*self.tica_concatenated[:, :2].T, self.flux.committor[self.dtrajs_concatenated], cmap='brg', ax=ax, mask=True, cbar_label=r'committor $\mathcal{S}_%d \to \mathcal{S}_%d$' % (start + 1, final + 1))
        fig.tight_layout()
        fig.savefig("onto_TICA.png")

