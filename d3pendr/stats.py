import numpy as np
from scipy import stats


def wasserstein_distance_and_direction(u_values, v_values, log10_transform):

    # modified from scipy.stats._cdf_distance

    u_sorter = np.argsort(u_values)
    v_sorter = np.argsort(v_values)

    all_values = np.concatenate((u_values, v_values))
    all_values.sort(kind='mergesort')

    # Compute the differences between pairs of successive values of u and v.
    deltas = np.diff(all_values)
    if log10_transform:
        deltas = np.log10(deltas + 1)

    # Get the respective positions of the values of u and v among the values of
    # both distributions.
    u_cdf_indices = u_values[u_sorter].searchsorted(all_values[:-1], 'right')
    v_cdf_indices = v_values[v_sorter].searchsorted(all_values[:-1], 'right')

    # Calculate the CDFs of u and v
    u_cdf = u_cdf_indices / u_values.size
    v_cdf = v_cdf_indices / v_values.size

    # Compute the value of the integral based on the CDFs.
    mag = np.multiply(u_cdf - v_cdf, deltas)

    # sum of abs is true wasserstein distance, sum accounts for direction of movement
    return np.sum(np.abs(mag)), np.sum(mag)


def wasserstein_test(u_values, v_values, bootstraps=999):
    # permutation test of wasserstein distance
    # based on the one outlined in https://github.com/cdowd/twosamples
    obs = stats.wasserstein_distance(u_values, v_values)

    # under null hypothesis the samples are drawn from the same distribution
    # so we can make expected wasserstein values by permuting values between
    # the two samples
    pool = np.concatenate([u_values, v_values])
    n = len(u_values)
    exp = []
    for _ in range(bootstraps):
        np.random.shuffle(pool)
        exp.append(stats.wasserstein_distance(pool[:n], pool[n:]))
    exp = np.array(exp)

    # bootstrap p value with pseudocount
    p = ((exp >= obs).sum() + 1) / (bootstraps + 1)
    return p


def tpe_stats(tpe_dist_cntrl, tpe_dist_treat, bootstraps=999, log10_transform=False):
    wass_dist, wass_dir = wasserstein_distance_and_direction(
        tpe_dist_cntrl, tpe_dist_treat, log10_transform
    )
    wass_pval = wasserstein_test(tpe_dist_cntrl, tpe_dist_treat, bootstraps)
    ks_stat, ks_pval = stats.ks_2samp(tpe_dist_cntrl, tpe_dist_treat)

    # combine the two tests using the harmonic mean of the p values
    try:
        hm_pval = stats.hmean([wass_pval, ks_pval])
    except ValueError:
        # if result is so sig that ks_2samp reports a pval of 0.0, hmean throws error
        # use the smallest number that can be represented as float64 as placeholder
        hm_pval = stats.hmean([wass_pval, np.finfo(np.float64).tiny])
    return wass_dist, wass_dir, wass_pval, ks_stat, ks_pval, hm_pval