from joblib import Parallel, delayed

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .io import MultiBamParser, bed12_iterator, write_output_bed
from .invs import intersect_spliced_invs
from .stats import tpe_stats


# columns in the output bed format
RESULTS_COLUMNS = [
    'chrom', 'start', 'end', 'gene_id', 'score', 'strand',
    'wass_dist', 'wass_dir', 'wass_pval', 'wass_fdr',
    'ks_stat',  'ks_pval', 'ks_fdr',
    'hm_pval', 'hm_fdr',
    'nreads_cntrl', 'nreads_treat'
]


def relative_tpe(gene_start, gene_end, aln_start, aln_end, strand):
    # make tpe position relative to gene start rather than 
    # genomic coordinates
    if strand == '+':
        return aln_end - gene_start
    else:
        return gene_end - aln_start


def get_tpe_distribution(mbam, chrom, gene_start, gene_end, strand,
                         gene_invs, min_read_overlap):
    tpe_distrib = []
    # fetch parsed alignments from bam files, filtering by strand
    for aln_start, aln_end, _, aln_invs in mbam.fetch(
            chrom, gene_start, gene_end, strand=strand):

        # calculate the fraction of the read alignment overlapping the
        # annotated gene
        aln_len = sum([e - s for s, e in aln_invs])
        i = intersect_spliced_invs(aln_invs, gene_invs)

        # only use the read if it overlaps with the reference annotation
        # by at least min_read_overlap fraction of its aligned length
        if i / aln_len >= min_read_overlap:
            tpe = relative_tpe(
                gene_start, gene_end,
                aln_start, aln_end,
                strand
            )
            tpe_distrib.append(tpe)

    return np.array(tpe_distrib)


def process_bed_records(bed_records, treat_bam_fns, cntrl_bam_fns,
                        min_read_overlap, bootstraps, log10_transform):
    results = []

    with MultiBamParser(cntrl_bam_fns) as cntrl_bam, \
            MultiBamParser(treat_bam_fns) as treat_bam:

        for chrom, start, end, gene_id, strand, invs in bed_records:
            cntrl_distrib = get_tpe_distribution(
                cntrl_bam, chrom, start, end, strand, invs,
                min_read_overlap
            )
            treat_distrib = get_tpe_distribution(
                treat_bam, chrom, start, end, strand, invs,
                min_read_overlap
            )
            nreads_cntrl = len(cntrl_distrib)
            nreads_treat = len(treat_distrib)

            if nreads_cntrl and nreads_treat:
                wass_dist, wass_dir, wass_pval, ks_stat, ks_pval, hm_pval = tpe_stats(
                    cntrl_distrib, treat_distrib,
                    bootstraps=bootstraps, log10_transform=log10_transform,
                )
                results.append([
                    chrom, start, end, gene_id, round(wass_dist), strand,
                    wass_dist, wass_dir, wass_pval, 1, # placeholder for wasserstein test fdr
                    ks_stat, ks_pval, 1, # placeholder for KS test fdr
                    hm_pval, 1, # placeholder for harmonic mean fdr
                    nreads_cntrl, nreads_treat
                ])

    results = pd.DataFrame(results, columns=RESULTS_COLUMNS)
    return results


def chunk_bed_records(bed_records, processes):
    # read the whole bed file
    bed_records = list(bed_records)
    nrecords = len(bed_records)
    n, r = divmod(nrecords, processes)
    split_points = ([0] + r * [n + 1] + (processes - r) * [n])
    split_points = np.cumsum(split_points)
    for i in range(processes):
        start = split_points[i]
        end = split_points[i + 1]
        yield bed_records[start: end]


def ref_guided_diff_tpe(bed_fn, treat_bam_fns, cntrl_bam_fns,
                        min_read_overlap, bootstraps,
                        log10_transform, processes):
    bed_it = bed12_iterator(bed_fn)
    if processes == 1:
        # run on main process
        results = process_bed_records(
            bed_it,
            treat_bam_fns,
            cntrl_bam_fns,
            min_read_overlap,
            bootstraps,
            log10_transform,log10_transform
        )
    else:
        args = (
            treat_bam_fns, cntrl_bam_fns, min_read_overlap,
            bootstraps, log10_transform
        )
        results = Parallel(n_jobs=processes)(
            delayed(process_bed_records)(bed_chunk, *args)
            for bed_chunk in chunk_bed_records(bed_it, processes)
        )
        results = pd.concat(results)
    _, results['wass_fdr'], *_ = multipletests(results.wass_pval, method='fdr_bh')
    _, results['ks_fdr'], *_ = multipletests(results.ks_pval, method='fdr_bh')
    _, results['hm_fdr'], *_ = multipletests(results.hm_pval, method='fdr_bh')
    return results