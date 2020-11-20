from bisect import bisect_right
from collections import defaultdict, Counter

import numpy as np
from scipy import ndimage as ndi
from scipy import signal


def relative_tpe(gene_start, gene_end, aln_start, aln_end, strand, read_end):
    # make tpe position relative to gene start rather than 
    # genomic coordinates
    if read_end == '3':
        if strand == '+':
            return aln_end - gene_start
        else:
            return gene_end - aln_start
    elif read_end == '5':
        if strand == '+':
            return aln_start - gene_start
        else:
            return gene_end - aln_end


def intersect(inv_a, inv_b):
    a_start, a_end = inv_a
    b_start, b_end = inv_b
    if a_end < b_start or a_start > b_end:
        return 0
    else:
        s = max(a_start, b_start)
        e = min(a_end, b_end)
        return e - s


def intersect_spliced_invs(invs_a, invs_b):
    score = 0
    invs_a = iter(invs_a)
    invs_b = iter(invs_b)
    a_start, a_end = next(invs_a)
    b_start, b_end = next(invs_b)
    while True:
        if a_end < b_start:
            try:
                a_start, a_end = next(invs_a)
            except StopIteration:
                break
        elif a_start > b_end:
            try:
                b_start, b_end = next(invs_b)
            except StopIteration:
                break
        else:
            score += intersect([a_start, a_end], [b_start, b_end])
            if a_end > b_end:
                try:
                    b_start, b_end = next(invs_b)
                except StopIteration:
                    break
            else:
                try:
                    a_start, a_end = next(invs_a)
                except StopIteration:
                    break
    return score


def get_tpe_distribution(mbam, chrom, gene_start, gene_end, strand,
                         min_read_overlap, read_strand, read_end,
                         paired_end_read, is_bam, as_coverage=False):
    tpe_distribs = []
    nreads = []

    # fetch parsed alignments from bam files, filtering by strand
    if is_bam:

        for sample in mbam.fetch(
                chrom, gene_start, gene_end,
                strand=strand, # the strand that the gene feature is on
                read_strand=read_strand, # the relative orientation of the reads to fetch
                pairs=paired_end_read): # which read end to fetch
            # sample is an iterator of parsed bam file alignments
            sample_distrib = []
            for aln_start, aln_end, aln_strand, aln_invs in sample:
                # calculate the fraction of the read alignment overlapping the
                # annotated gene
                aln_len = sum([e - s for s, e in aln_invs])
                i = intersect_spliced_invs(aln_invs, [(gene_start, gene_end),])

                # only use the read if it overlaps with the reference annotation
                # by at least min_read_overlap fraction of its aligned length
                if i / aln_len >= min_read_overlap:
                    tpe = relative_tpe(
                        gene_start, gene_end,
                        aln_start, aln_end,
                        strand, read_end
                    )
                    sample_distrib.append(tpe)
            if as_coverage:
                ln = gene_end - gene_start
                sample_distrib = np.bincount(sample_distrib, minlength=ln)[:ln]
                n = sum(sample_distrib)
            else:
                n = len(sample_distrib)
            tpe_distribs.append(np.array(sample_distrib))
            nreads.append(n)

    else:
        for sample in mbam.coverage(
                chrom, gene_start, gene_end,
                strand=fetch_strand, normalise=False):
            # sample must be count data
            if np.mod(sample, 1).any():
                raise ValueError('Bigwig values are not discrete')
            sample = sample.astype(np.int)
            # sample is a hist of coverage over the gene from the bigwig file
            # it needs inverting for negative strand genes
            if strand == '-':
                sample = sample[::-1]
            if not as_coverage:
                sample_distrib = np.repeat(np.arange(gene_end - gene_start), sample)
                n = len(sample_distrib)
            else:
                n = sum(sample_distrib)
            tpe_distribs.append(sample_distrib)
            nreads.append(n)

    nreads = np.array(nreads)
    return tpe_distribs, nreads


def find_tpe_peaks_in_coverage_bincounts(tpe_cov, sigma, min_exprs, min_dist):
    tpe_cov = ndi.convolve(
        tpe_cov, np.ones(sigma), mode='constant', cval=0
    )
    # in case there is a flat peak at the start/end of the cov
    tpe_cov = np.insert(tpe_cov, [0, len(tpe_cov)], 0, 0)

    tpes, _ = signal.find_peaks(
        tpe_cov, height=min_exprs, distance=min_dist
    )
    return tpes - 1


def argrelmin_left_on_flat(arr, order):
    idx = []
    for i in range(order, len(arr) - order):
        if np.all(arr[i] < arr[i - order: i]) and np.all(arr[i] <= arr[i + 1: i + 1 + order]):
            idx.append(i)
    return np.array(idx)


def cluster_by_endpoint(endpoints, conds, sigma):
    offset = endpoints.min() - sigma * 3
    endpoints_scaled = endpoints - offset
    endpoints_max = endpoints_scaled.max()
    endpoints_dist = np.bincount(
        endpoints_scaled,
        minlength=endpoints_max + sigma * 3
    ).astype('f')
    endpoints_dist = ndi.filters.gaussian_filter1d(
        endpoints_dist, sigma=sigma, mode='constant', cval=0
    )

    # find local minima in three prime positions
    cut_idx = argrelmin_left_on_flat(endpoints_dist, order=sigma)
    cut_idx = cut_idx + offset

    # classify alignments in relation to local minima
    cluster_idx = defaultdict(list)
    cluster_cond_count = defaultdict(list)
    for pos, c in zip(endpoints, conds):
        i = bisect_right(cut_idx, pos)
        cluster_idx[i].append(pos)
        cluster_cond_count[i].append(c)

    # get actual start/end points of cluster and cluster counts for each cond
    clusters = {}
    for i, pos in cluster_idx.items():
        inv = (min(pos), max(pos))
        clusters[inv] = Counter(cluster_cond_count[i])
    return clusters


def get_apa_tpes(cntrl_distrib, nreads_cntrl,
                 treat_distrib, nreads_treat,
                 sigma, min_count, min_rel_change):
    endpoints = np.concatenate([
        *cntrl_distrib, *treat_distrib
    ])
    conds = np.repeat(
        [0, 1],
        [nreads_cntrl, nreads_treat]
    )
    assert len(endpoints) == len(conds)
    tpes = cluster_by_endpoint(endpoints, conds, sigma)

    for (start, end), counts in tpes.items():
        cntrl_count = counts[0]
        treat_count = counts[1]
        if (cntrl_count + treat_count) >= min_count:
            cntrl_frac = cntrl_count / nreads_cntrl
            treat_frac = treat_count / nreads_treat
            relative_change = treat_frac - cntrl_frac
            if np.abs(relative_change) >= min_rel_change:
                yield (start, end, cntrl_count, treat_count,
                       cntrl_frac, treat_frac, relative_change)