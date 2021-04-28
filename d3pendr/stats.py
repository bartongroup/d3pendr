import itertools as it
import numpy as np
from scipy import stats


def wasserstein_distance_and_direction(u_values, v_values):

    # modified from scipy.stats._cdf_distance

    u_sorter = np.argsort(u_values)
    v_sorter = np.argsort(v_values)

    all_values = np.concatenate((u_values, v_values))
    all_values = np.sort(all_values)

    # Compute the differences between pairs of successive values of u and v.
    deltas = np.diff(all_values)

    # Get the respective positions of the values of u and v among the values of
    # both distributions.
    u_cdf_indices = np.searchsorted(u_values[u_sorter], all_values[:-1], 'right')
    v_cdf_indices = np.searchsorted(v_values[v_sorter], all_values[:-1], 'right')

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


def wasserstein_silhouette_test(cntrl, treat, bootstraps=999,
                                fit_skewnorm=True, fit_gamma=True,
                                random_state=None):
    rs = np.random.RandomState(random_state)
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
    n = split_idx[len(cntrl) - 1]
    pool = np.concatenate(pool)
    exp_wass = []
    exp_sil = []
    for _ in range(bootstraps):
        rs.shuffle(pool)
        exp_wass.append(stats.wasserstein_distance(pool[:n], pool[n:]))
        shuf_samps = np.array_split(pool, split_idx)
        exp_dist_mat = pairwise_wasserstein(shuf_samps)
        exp_sil.append(silhouette_score(exp_dist_mat, labels))
    exp_wass = np.array(exp_wass)
    exp_sil = np.array(exp_sil)

    if fit_gamma:
        g_fit = stats.gamma.fit(exp_wass)
        wass_pval = stats.gamma.sf(wass_dist, *g_fit)
    else:
        wass_pval = ((exp >= wass_dist).sum() + 1) / (bootstraps + 1)

    if fit_skewnorm:
        s_fit = stats.skewnorm.fit(exp_sil)
        sil_pval = stats.skewnorm.sf(sil, *s_fit)
    else:
        sil_pval = ((sil <= exp_sil).sum() + 1) / (bootstraps + 1)
    return wass_dist, wass_dir, wass_pval, sil, sil_pval


def wasserstein_test(u_values, v_values, bootstraps=999,
                     fit_gamma=True, random_state=None):

    rs = np.random.RandomState(random_state)
    # permutation test of wasserstein distance
    # based on the one outlined in https://github.com/cdowd/twosamples
    wass_dist, wass_dir = wasserstein_distance_and_direction(u_values, v_values)

    # under null hypothesis the samples are drawn from the same distribution
    # so we can make expected wasserstein values by permuting values between
    # the two samples
    pool = np.concatenate([u_values, v_values])
    n = len(u_values)
    exp = []
    for _ in range(bootstraps):
        rs.shuffle(pool)
        exp.append(stats.wasserstein_distance(pool[:n], pool[n:]))
    exp = np.array(exp)

    if fit_gamma:
        # fit a gamma distribution to the expected distances
        g_fit = stats.gamma.fit(exp_wass)
        # compute p value using survival function
        wass_p_val = stats.gamma.sf(wass_dist, *g_fit)
    else:
        # bootstrap p value with pseudocount
        p_val = ((exp >= wass_dist).sum() + 1) / (bootstraps + 1)
    return wass_dist, wass_dir, p_val


def d3pendr_stats(tpe_dist_cntrl, tpe_dist_treat,
                  bootstraps=999, fit_gamma=True,
                  sil_test=True, fit_skewnorm=True,
                  random_state=None):
    if len(tpe_dist_cntrl) == 1 or not sil_test:
        wass_dist, wass_dir, wass_pval = wasserstein_test(
            np.concatenate(tpe_dist_cntrl),
            np.concatenate(tpe_dist_treat),
            bootstraps, fit_gamma, random_state
        )
        sil, sil_pval = np.nan, np.nan
    else:
        wass_dist, wass_dir, wass_pval, sil, sil_pval = wasserstein_silhouette_test(
            tpe_dist_cntrl, tpe_dist_treat, bootstraps,
            fit_gamma, fit_skewnorm, random_state
        )

    return wass_dist, wass_dir, wass_pval, sil, sil_pval