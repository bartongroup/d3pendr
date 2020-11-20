import click

from .output import (
    write_wass_test_output_bed,
    write_apa_site_bed,
)
from .d3pendr import run_d3pendr


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
@click.option('--min-reads-per-rep', default=5)
@click.option('--extend-gene-five-prime', default=0)
@click.option('--use-5utr/--ignore-5utr', default=True)
@click.option('--extend-gene-three-prime', default=0)
@click.option('--use-locus-tag/--use-gene-id-tag', default=False)
@click.option('--max-terminal-intron-size', default=100_000)
@click.option('--min-terminal-exon-size', default=30)
@click.option('--bootstraps', default=999)
@click.option('--threshold', default=0.05)
@click.option('--use-gamma-model/--no-model', default=True)
@click.option('--test-homogeneity/--no-test-homogeneity', default=False)
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
            use_locus_tag,
            max_terminal_intron_size, min_terminal_exon_size,
            bootstraps, threshold, use_gamma_model, test_homogeneity,
            tpe_cluster_sigma, min_tpe_reads, min_tpe_fractional_change,
            processes):
    '''
    d3pendr: Differential 3' End analysis of Nanopore Direct RNAs

    Identifies differential polyadenylation events in either an
    annotation dependent manner, using
    a permutation test of Wasserstein distance between 3' end 
    distributions and a KS test.

    Outputs bed6 format with extra columns.
    '''

    results, apa_sites = run_d3pendr(
        annotation_gtf_fn,
        treatment_fns, control_fns,
        read_strand, read_end, paired_end_read,
        min_read_overlap, min_reads_per_rep,
        extend_gene_five_prime, use_5utr,
        extend_gene_three_prime,
        use_locus_tag,
        max_terminal_intron_size, min_terminal_exon_size,
        bootstraps, threshold,
        use_gamma_model, test_homogeneity,
        write_apa_sites,
        tpe_cluster_sigma, min_tpe_reads, min_tpe_fractional_change,
        processes
    )
    output_bed = f'{output_prefix}.apa_results.bed'
    write_wass_test_output_bed(output_bed, results)
    if write_apa_sites:
        apa_site_bed = f'{output_prefix}.apa_sites.bed'
        write_apa_site_bed(apa_site_bed, apa_sites)


if __name__ == '__main__':
    d3pendr()