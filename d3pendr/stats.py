import itertools as it
import numpy as np
from scipy import stats


def wasserstein_distance_and_direction(u_values, v_values):

    # modified from scipy.stats._cdf_distance

    u_sorter = np.argsort(u_values)
    v_sorter = np.argsort(v_values)

    all_values = np.concatenate((u_values, v_values))
    all_values.sort(kind='mergesort')

    # Compute the differences between pairs of successive values of u and v.
    deltas = np.diff(all_values)

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
    wass_dist, wass_dir = wasserstein_distance_and_direction(u_values, v_values)

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
    p_val = ((exp >= wass_dist).sum() + 1) / (bootstraps + 1)
    return wass_dist, wass_dir, p_val


def wasserstein_test_with_replicates(u_values, v_values,
                                     bootstraps=999, threshold=0.05):
    # first pool replicates and perform normal wass test
    u_pooled = np.concatenate(u_values)
    v_pooled = np.concatenate(v_values)
    wass_dist, wass_dir, p_val = wasserstein_test(u_pooled, v_pooled, bootstraps=bootstraps)

    # we only bother testing for homogeneity of replicates if the test between
    # replicates is significant
    if p_val <= threshold:

        # produce pairwise distances between replicates for both conditions
        u_self_wass = np.array(
            [stats.wasserstein_distance(u_i, u_j)
             for u_i, u_j in it.combinations(u_values, r=2)]
        )
        v_self_wass = np.array(
            [stats.wasserstein_distance(v_i, v_j)
             for v_i, v_j in it.combinations(v_values, r=2)]
        )

        # samples are considered homogeneous if the distance between the two
        # conds is greater than all the pairwise distances within conds.
        is_homogeneous = 1 - int((u_self_wass >= wass_dist).any() |
                                 (v_self_wass >= wass_dist).any())
        if not is_homogeneous:
            p_val = 1
    return wass_dist, wass_dir, p_val


def tpe_stats(tpe_dist_cntrl, tpe_dist_treat, bootstraps=999, threshold=0.05):
    if len(tpe_dist_cntrl) == 1:
        wass_dist, wass_dir, wass_pval = wasserstein_test(
            tpe_dist_cntrl[0], tpe_dist_treat[0], bootstraps
        )
    else:
        wass_dist, wass_dir, wass_pval = wasserstein_test_with_replicates(
            tpe_dist_cntrl, tpe_dist_treat, bootstraps, threshold
        )

    return wass_dist, wass_dir, wass_pval