import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed

from .multibam import MultiBamParser, MultiBigWigParser, bam_or_bw
from .gtf import gtf_iterator, chunk_gtf_records
from .tpe import get_tpe_distribution, get_apa_tpes
from .stats import d3pendr_stats


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
        

def _process_gtf_records(gtf_records, treat_bam_fns, cntrl_bam_fns,
                         read_strand, read_end, paired_end_read,
                         min_read_overlap, min_reads,
                         bootstraps, threshold,
                         use_gamma_model, test_homogeneity,
                         is_bam, find_apa_tpe_sites, tpe_sigma,
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
                wass_dist, wass_dir, wass_pval = d3pendr_stats(
                    cntrl_distribs, treat_distribs,
                    bootstraps=bootstraps, threshold=threshold,
                    use_gamma_model=use_gamma_model,
                    test_homogeneity=test_homogeneity,
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


def run_d3pendr(gtf_fn, treat_bam_fns, cntrl_bam_fns,
                read_strand, read_end, paired_end_read,
                min_read_overlap, min_reads_per_cond,
                extend_gene_five_prime, use_5utr,
                extend_gene_three_prime,
                by_locus,
                max_terminal_intron_size,
                min_terminal_exon_size,
                bootstraps, threshold,
                use_gamma_model, test_homogeneity,
                find_tpe_sites, tpe_sigma,
                tpe_min_reads, tpe_min_rel_change,
                processes):

    filetype = bam_or_bw(*treat_bam_fns, *cntrl_bam_fns)

    args = (
        treat_bam_fns, cntrl_bam_fns,
        read_strand, read_end, paired_end_read,
        min_read_overlap, min_reads_per_cond,
        bootstraps, threshold,
        use_gamma_model, test_homogeneity,
        filetype, find_tpe_sites, tpe_sigma,
        tpe_min_reads, tpe_min_rel_change,
    )
    gtf_it = gtf_iterator(
        gtf_fn, extend_gene_five_prime, use_5utr, extend_gene_three_prime,
        by_locus, max_terminal_intron_size, min_terminal_exon_size,
    )
    if processes == 1:
        # run on main process
        results, tpa_apa_results = _process_gtf_records(
            gtf_it, *args
        )
    else:
        results = Parallel(n_jobs=processes)(
            delayed(_process_gtf_records)(gtf_chunk, *args)
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