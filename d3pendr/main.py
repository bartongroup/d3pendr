import click

from .io import write_output_bed, write_apa_site_bed
from .ref_guided import ref_guided_diff_tpe


@click.command()
@click.option('-t', '--treatment-fns', required=True, multiple=True)
@click.option('-c', '--control-fns', required=True, multiple=True)
@click.option('-o', '--output-prefix', required=True)
@click.option('--write-apa-sites/--no-apa-sites', required=False, default=True)
@click.option('-a', '--annotation-gtf-fn', required=True)
@click.option('--read-strand', type=click.Choice(['same', 'opposite', 'unstranded']), default='same')
@click.option('--read-end', type=click.Choice(['3', '5']), default='3')
@click.option('--paired-end-read', type=click.Choice(['1', '2', 'both', 'single']), default='single')
@click.option('--min-read-overlap', default=0.2)
@click.option('--min-reads-per-rep', default=10)
@click.option('--extend-gene-five-prime', default=0)
@click.option('--use-5utr/--ignore-5utr', default=True)
@click.option('--extend-gene-three-prime', default=0)
@click.option('--bootstraps', default=1999)
@click.option('--threshold', default=0.05)
@click.option('--tpe-cluster-sigma', default=15)
@click.option('--min-tpe-reads', default=5)
@click.option('--min-tpe-fractional-change', default=0.1)
@click.option('-p', '--processes', default=4)
def d3pendr(treatment_fns, control_fns,
            output_prefix, write_apa_sites,
            annotation_gtf_fn,
            read_strand, read_end, paired_end_read,
            min_read_overlap, min_reads_per_rep,
            extend_gene_five_prime, use_5utr,
            extend_gene_three_prime,
            bootstraps, threshold, tpe_cluster_sigma,
            min_tpe_reads, min_tpe_fractional_change,
            processes):
    '''
    d3pendr: Differential 3' End analysis of Nanopore Direct RNAseq

    Identifies differential polyadenylation events in either an
    annotation dependent manner, using
    a permutation test of Wasserstein distance between 3' end 
    distributions and a KS test.

    Outputs bed6 format with extra columns.
    '''

    if paired_end_read == 'single':
        # for options single or both, all reads are used
        paired_end_read = 'both'

    results, apa_sites = ref_guided_diff_tpe(
        annotation_gtf_fn,
        treatment_fns, control_fns,
        read_strand, read_end, paired_end_read,
        min_read_overlap, min_reads_per_rep,
        extend_gene_five_prime, use_5utr,
        extend_gene_three_prime,
        bootstraps, threshold, write_apa_sites,
        tpe_cluster_sigma, min_tpe_reads, min_tpe_fractional_change,
        processes
    )
    output_bed = f'{output_prefix}.apa_results.bed'
    write_output_bed(output_bed, results)
    if write_apa_sites:
        apa_site_bed = f'{output_prefix}.apa_sites.bed'
        write_apa_site_bed(apa_site_bed, apa_sites)


if __name__ == '__main__':
    d3pendr()