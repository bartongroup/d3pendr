import itertools as it
from bisect import bisect_right
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed

from .io import (
    MultiBamParser, MultiBigWigParser, bam_or_bw,
    gtf_iterator, write_output_bed
)
from .invs import relative_tpe, intersect_spliced_invs
from .stats import tpe_stats


# columns in the output bed format
RESULTS_COLUMNS = [
    'chrom', 'start', 'end', 'gene_id', 'score', 'strand',
    'wass_dist', 'wass_dir', 'wass_pval', 'wass_fdr',
    'nreads_cntrl', 'nreads_treat'
]

TPE_APA_RESULTS_COLUMNS = [
    'chrom', 'start', 'end', 'gene_id', 'strand',
    'nreads_cntrl', 'nreads_treat',
    'frac_cntrl', 'frac_treat', 'relative_change'
]


def get_tpe_distribution(mbam, chrom, gene_start, gene_end, strand,
                         min_read_overlap, read_strand, read_end,
                         paired_end_read, is_bam):
    tpe_distribs = []
    nreads = []

    if read_strand == 'opposite':
        fetch_strand = '+' if strand == '-' else '-'
    elif read_strand == 'unstranded':
        fetch_strand = None
    else:
        fetch_strand = strand

    # fetch parsed alignments from bam files, filtering by strand
    if is_bam:

        for sample in mbam.fetch(
                chrom, gene_start, gene_end,
                strand=fetch_strand, pairs=paired_end_read):
            # sample is an iterator of parsed bam file alignments
            
            sample_distrib = []
            n = 0
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
                    n += 1
            tpe_distribs.append(np.array(sample_distrib))
            nreads.append(n)

    else:
        for sample in mbam.fetch(
                chrom, gene_start, gene_end, strand=fetch_strand):
            # sample must be count data
            sample[np.isnan(sample)] = 0
            if np.mod(sample, 1).any():
                raise ValueError('Bigwig values are not discrete')
            sample = sample.astype(np.int)
            # sample is a hist of coverage over the gene from the bigwig file
            # it needs inverting for negative strand genes
            if strand == '-':
                sample = sample[::-1]
            sample_distrib = np.repeat(np.arange(gene_end - gene_start), sample)
            tpe_distribs.append(sample_distrib)
            nreads.append(len(sample_distrib))

    nreads = np.array(nreads)
    return tpe_distribs, nreads


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
    endpoints_dist = gaussian_filter1d(
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
        

def process_gtf_records(gtf_records, treat_bam_fns, cntrl_bam_fns,
                        read_strand, read_end, paired_end_read,
                        min_read_overlap, min_reads,
                        bootstraps, threshold, is_bam,
                        find_apa_tpe_sites, tpe_sigma,
                        tpe_min_reads, tpe_min_rel_change):
    results = []
    tpe_apa_res = []

    parser = MultiBamParser if is_bam else MultiBigWigParser

    with parser(cntrl_bam_fns) as cntrl_bam, parser(treat_bam_fns) as treat_bam:

        for chrom, start, end, gene_id, strand in gtf_records:
            cntrl_distribs, nreads_cntrl = get_tpe_distribution(
                cntrl_bam, chrom, start, end, strand,
                min_read_overlap, read_strand, read_end,
                paired_end_read, is_bam
            )
            treat_distribs, nreads_treat = get_tpe_distribution(
                treat_bam, chrom, start, end, strand,
                min_read_overlap, read_strand, read_end,
                paired_end_read, is_bam
            )
                    

            if (nreads_cntrl >= min_reads).all() and (nreads_treat >= min_reads).all():
                wass_dist, wass_dir, wass_pval = tpe_stats(
                    cntrl_distribs, treat_distribs,
                    bootstraps=bootstraps, threshold=threshold,
                )
                results.append([
                    chrom, start, end, gene_id, round(wass_dist), strand,
                    wass_dist, wass_dir, wass_pval, 1, # placeholder for wasserstein test fdr
                    nreads_cntrl, nreads_treat
                ])

                if find_apa_tpe_sites and wass_pval <= threshold:
                    tpes = get_apa_tpes(
                        cntrl_distribs, sum(nreads_cntrl),
                        treat_distribs, sum(nreads_treat),
                        tpe_sigma, tpe_min_reads, tpe_min_rel_change
                    )
                    for t in tpes:
                        tpe_start, tpe_end, *res = t
                        # revert tpe coords to absolute
                        if strand == '+':
                            tpe_start = tpe_start + start
                            tpe_end = tpe_end + start
                        elif strand == '-':
                            tpe_start, tpe_end = tpe_end, tpe_start
                            tpe_start = end - tpe_start
                            tpe_end = end - tpe_end

                        tpe_apa_res.append([
                            chrom, tpe_start, tpe_end, gene_id, strand, *res
                        ])

    results = pd.DataFrame(results, columns=RESULTS_COLUMNS)

    if find_apa_tpe_sites:
        tpe_apa_res = pd.DataFrame(tpe_apa_res, columns=TPE_APA_RESULTS_COLUMNS)
    else:
        tpe_apa_res = None
    return results, tpe_apa_res


def chunk_gtf_records(gtf_records, processes):
    # read the whole gtf file
    gtf_records = list(gtf_records)
    nrecords = len(gtf_records)
    n, r = divmod(nrecords, processes)
    split_points = ([0] + r * [n + 1] + (processes - r) * [n])
    split_points = np.cumsum(split_points)
    for i in range(processes):
        start = split_points[i]
        end = split_points[i + 1]
        yield gtf_records[start: end]


def ref_guided_diff_tpe(gtf_fn, treat_bam_fns, cntrl_bam_fns,
                        read_strand, read_end, paired_end_read,
                        min_read_overlap, min_reads_per_cond,
                        extend_gene_five_prime, use_5utr,
                        extend_gene_three_prime,
                        bootstraps, threshold,
                        find_tpe_sites, tpe_sigma,
                        tpe_min_reads, tpe_min_rel_change,
                        processes):

    filetypes = [bam_or_bw(fn) for fn in it.chain(treat_bam_fns, cntrl_bam_fns)]
    assert all([t == filetypes[0] for t in filetypes])
    filetype = filetypes[0]
    args = (
        treat_bam_fns, cntrl_bam_fns,
        read_strand, read_end, paired_end_read,
        min_read_overlap, min_reads_per_cond,
        bootstraps, threshold, filetype,
        find_tpe_sites, tpe_sigma,
        tpe_min_reads, tpe_min_rel_change,
    )
    gtf_it = gtf_iterator(
        gtf_fn, extend_gene_five_prime, use_5utr, extend_gene_three_prime
    )
    if processes == 1:
        # run on main process
        results, tpa_apa_results = process_gtf_records(
            gtf_it, *args
        )
    else:
        results = Parallel(n_jobs=processes)(
            delayed(process_gtf_records)(gtf_chunk, *args)
            for gtf_chunk in chunk_gtf_records(gtf_it, processes)
        )
        results, tpe_apa_results = zip(*results)
        results = pd.concat(results)
        if find_tpe_sites:
            tpe_apa_results = pd.concat(tpe_apa_results)
        else:
            tpe_apa_results = None

    _, results['wass_fdr'], *_ = multipletests(results.wass_pval, method='fdr_bh')

    if find_tpe_sites:
        tpe_apa_results = tpe_apa_results[
            tpe_apa_results.gene_id.isin(results.query(f'wass_fdr <= {threshold}').gene_id)
        ]
    return results, tpe_apa_results