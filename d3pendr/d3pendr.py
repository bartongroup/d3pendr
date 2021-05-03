import dataclasses
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from joblib import Parallel, delayed

from .multibam import MultiBamParser, MultiBigWigParser, bam_or_bw
from .gtf import GTFGeneBoundaryIterator, chunk_gtf_records, GTFRecord
from .tpe import get_tpe_distribution
from .stats import d3pendr_stats, WassTestResults


@dataclasses.dataclass
class D3pendrResults(WassTestResults, GTFRecord):
    '''class for handling all d3pendr results for a record'''

    @staticmethod
    def from_records(gtf_record, wasstest_results):
        return D3pendrResults(
            **dataclasses.asdict(gtf_record),
            **dataclasses.asdict(wasstest_results)
        )


def _process_gtf_records(gtf_records, bam_opts, stat_opts, random_seed, is_bam):

    random_state = np.random.default_rng(random_seed)
    results = []
    tpe_apa_res = []

    parser = MultiBamParser if is_bam else MultiBigWigParser

    with parser(bam_opts.control_fns) as cntrl_bam, parser(bam_opts.treatment_fns) as treat_bam:

        for record in gtf_records:
            cntrl_distribs, nreads_cntrl = get_tpe_distribution(
                cntrl_bam,
                record.chrom, record.start, record.end, record.strand,
                bam_opts.min_read_overlap,
                bam_opts.read_strand, bam_opts.read_end,
                bam_opts.paired_end_read, is_bam
            )
            treat_distribs, nreads_treat = get_tpe_distribution(
                treat_bam,
                record.chrom, record.start, record.end, record.strand,
                bam_opts.min_read_overlap,
                bam_opts.read_strand, bam_opts.read_end,
                bam_opts.paired_end_read, is_bam
            )

            min_reads = bam_opts.min_reads_per_rep
            if (nreads_cntrl >= min_reads).all() and (nreads_treat >= min_reads).all():
                wass_res = d3pendr_stats(
                    cntrl_distribs, treat_distribs,
                    random_state=random_state,
                    **dataclasses.asdict(stat_opts)
                )
                results.append(D3pendrResults.from_records(record, wass_res))

    results = pd.DataFrame(results)
    return results


def run_d3pendr(opts):

    filetype = bam_or_bw(*opts.bam.treatment_fns, *opts.bam.control_fns)
    gtf_it = GTFGeneBoundaryIterator(**dataclasses.asdict(opts.gtf))
    if opts.processes == 1:
        # run on main process
        results, tpa_apa_results = _process_gtf_records(
            gtf_it, opts.bam, opts.stats, opts.random_seed, filetype
        )
    else:
        results = Parallel(n_jobs=opts.processes)(
            delayed(_process_gtf_records)(
                gtf_chunk, opts.bam, opts.stats, opts.random_seed, filetype
            ) for gtf_chunk in chunk_gtf_records(gtf_it, opts.processes)
        )
        results = pd.concat(results)

    _, results['wass_fdr'], *_ = multipletests(results.wass_pval, method='fdr_bh')
    _, results['sil_fdr'], *_ = multipletests(results.sil_pval, method='fdr_bh')

    return results