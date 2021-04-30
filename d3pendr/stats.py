import dataclasses
import numpy as np
from scipy import stats


@dataclasses.dataclass
class WassTestResults:
    '''Class for handling results of wasserstein_test and wasserstein_silhouette_test'''
    nreads_cntrl: int
    nreads_treat: int
    wass_dist: float
    wass_dir: float
    wass_pval: float
    wass_fdr: float = np.nan
    silhouette: float = np.nan
    sil_pval: float = np.nan
    sil_fdr: float = np.nan


def wasserstein_distance_and_direction(u_values, v_values):

    # modified from scipy.stats._cdf_distance

    u_values.sort()
    v_values.sort()

    all_values = np.concatenate((u_values, v_values))
    all_values.sort()

    # Compute the differences between pairs of successive values of u and v.
    deltas = np.diff(all_values)

    # Get the respective positions of the values of u and v among the values of
    # both distributions.
    u_cdf_indices = np.searchsorted(u_values, all_values[:-1], 'right')
    v_cdf_indices = np.searchsorted(v_values, all_values[:-1], 'right')

    # Calculate the CDFs of u and v
    u_cdf = u_cdf_indices / u_values.size
    v_cdf = v_cdf_indices / v_values.size

    # Compute the value of the integral based on the CDFs.
    mag = np.multiply(u_cdf - v_cdf, deltas)

    # sum of abs is true wasserstein distance, sum accounts for direction of movement
    return np.sum(np.abs(mag)), np.sum(mag)


def pairwise_wasserstein(samples):
    dist_mat = np.empty(shape=(len(samples), len(samples)))
    for i, j in zip(*np.triu_indices(len(samples))):
        if i == j:
            dist_mat[i, j] = 0
        else:
            wass, _ = wasserstein_distance_and_direction(samples[i], samples[j])
            dist_mat[i, j] = wass
            dist_mat[j, i] = wass
    return dist_mat


def silhouette_score(X, labels):
    assert np.array_equal(np.sort(labels), labels)
    # renumber from zero
    unique_labels = np.unique(labels, return_inverse=True)
    n_samples = len(labels)
    n_labels = len(unique_labels)

    idx = np.indices(X.shape).reshape(2, -1)
    sum_ = np.zeros((n_samples, n_labels))
    count = np.zeros((n_samples,  n_labels))
    for i, j, r, c, d in zip(*idx, *labels[idx], X.ravel()):
        if i != j:
            sum_[i, c] += d
            count[i, c] += 1
    mean = (sum_ / count).tolist()
    a = np.empty(n_samples)
    b = np.empty(n_samples)
    for i, (lab, vals) in enumerate(zip(labels, mean)):
        a[i] = vals.pop(lab)
        b[i] = min(vals)
    return np.mean((b - a) / np.maximum(a, b))


def wasserstein_silhouette_test(cntrl, treat, bootstraps=999,
                                fit_gamma=True, fit_skewnorm=True,
                                random_state=None):
    if not isinstance(random_state, np.random.Generator):
        random_state = np.random.default_rng(random_state)
    wass_dist, wass_dir = wasserstein_distance_and_direction(
        np.concatenate(cntrl),
        np.concatenate(treat)
    )
    labels = np.repeat([0, 1], [len(cntrl), len(treat)])
    pool = cntrl + treat
    obs_dist_mat = pairwise_wasserstein(pool)
    sil = silhouette_score(obs_dist_mat, labels)
    samp_lns = [len(samp) for samp in pool]
    split_idx = np.cumsum(samp_lns[:-1])
    nreads_cntrl = split_idx[len(cntrl) - 1]
    nreads_treat = sum(samp_lns) - nreads_cntrl
    pool = np.concatenate(pool)
    exp_wass = []
    exp_sil = []
    for _ in range(bootstraps):
        samp = random_state.choice(pool, size=len(pool), replace=True)
        exp_wass.append(wasserstein_distance_and_direction(
            samp[:nreads_cntrl], samp[nreads_cntrl:]
        )[0])
        shuf_samps = np.array_split(samp, split_idx)
        exp_dist_mat = pairwise_wasserstein(shuf_samps)
        exp_sil.append(silhouette_score(exp_dist_mat, labels))
    exp_wass = np.array(exp_wass)
    exp_sil = np.array(exp_sil)

    if fit_gamma:
        gfit = stats.gamma.fit(exp_wass)
        gamma = stats.gamma(*gfit)
        # compute p value using survival function
        wass_pval = gamma.sf(wass_dist)
    else:
        wass_pval = ((exp_wass >= wass_dist).sum() + 1) / (bootstraps + 1)

    if fit_skewnorm:
        sfit = stats.skewnorm.fit(exp_sil)
        skn = stats.skewnorm(*sfit)
        # compute p value using survival function
        sil_pval = skn.sf(sil)
    else:
        sil_pval = ((exp_sil >= sil).sum() + 1) / (bootstraps + 1)
        
    res = WassTestResults(
        nreads_cntrl, nreads_treat,
        wass_dist, wass_dir,
        wass_pval, 1, # placeholder for wass_fdr
        sil, sil_pval, 1, # placeholder
    )
    return res


def wasserstein_test(cntrl, treat, bootstraps=999,
                     fit_gamma=True, random_state=None):
    if not isinstance(random_state, np.random.Generator):
        random_state = np.random.default_rng(random_state)
    # permutation test of wasserstein distance
    # based on the one outlined in https://github.com/cdowd/twosamples
    wass_dist, wass_dir = wasserstein_distance_and_direction(cntrl, treat)

    # under null hypothesis the samples are drawn from the same distribution
    # so we can make expected wasserstein values by permuting values between
    # the two samples
    pool = np.concatenate([cntrl, treat])
    nreads_cntrl = len(cntrl)
    nreads_treat = len(treat)
    exp = []
    for _ in range(bootstraps):
        samp = random_state.choice(pool, size=len(pool), replace=True)
        exp.append(wasserstein_distance_and_direction(
            samp[:nreads_cntrl], samp[nreads_cntrl:]
        )[0])
    exp = np.array(exp)

    if fit_gamma:
        # fit a gamma distribution to the expected distances
        gfit = stats.gamma.fit(exp)
        gamma = stats.gamma(*gfit)
        # compute p value using survival function
        wass_pval = gamma.sf(wass_dist)
    else:
        # bootstrap p value with pseudocount
        wass_pval = ((exp >= wass_dist).sum() + 1) / (bootstraps + 1)

    res = WassTestResults(
        nreads_cntrl, nreads_treat,
        wass_dist, wass_dir, wass_pval, 1, # placeholder for wass_fdr
    )
    return res


def d3pendr_stats(tpe_dist_cntrl,
                  tpe_dist_treat,
                  bootstraps=999,
                  wass_fit_gamma=True,
                  silhouette_test=True,
                  sil_fit_skewnorm=True,
                  random_state=None,
                  **kwargs):
    if len(tpe_dist_cntrl) == 1 or not silhouette_test:
        wass_res = wasserstein_test(
            np.concatenate(tpe_dist_cntrl),
            np.concatenate(tpe_dist_treat),
            bootstraps, wass_fit_gamma, random_state
        )
    else:
        wass_res = wasserstein_silhouette_test(
            tpe_dist_cntrl, tpe_dist_treat, bootstraps,
            wass_fit_gamma, sil_fit_skewnorm, random_state
        )

    return wass_res