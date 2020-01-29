import click

from .io import write_output_bed
from .ref_guided import ref_guided_diff_tpe


def annot_free_diff_tpe(treat_bam_fns, cntrl_bam_fns,
                        min_read_overlap, bootstraps,
                        log10_transform, processes):
    raise NotImplementedError('TODO!')
    

@click.command()
@click.option('-t', '--treatment-bams', required=True, multiple=True)
@click.option('-c', '--control-bams', required=True, multiple=True)
@click.option('-o', '--output-bed', required=True)
@click.option('-a', '--annotation-bed12', required=False, default=None)
@click.option('--min-read-overlap', default=0.2)
@click.option('--bootstraps', default=999)
@click.option('--log10-transform', default=False)
@click.option('-p', '--processes', default=4)
def d3pendr(treatment_bams, control_bams, output_bed, annotation_bed12,
            min_read_overlap, bootstraps, log10_transform, processes):
    '''
    d3pendr: Differential 3' End analysis of Nanopore Direct RNAseq

    Identifies differential polyadenylation events in either an
    annotation dependent or (TODO) reference free manner, using
    a permutation test of Wasserstein distance between 3' end 
    distributions and a KS test.

    Outputs bed6 format with extra columns.
    '''
    if annotation_bed12 is None:
        results = annot_free_diff_tpe(
            treatment_bams, control_bams,
            min_read_overlap, bootstraps,
            log10_transform,
            processes
        )
    else:
        results = ref_guided_diff_tpe(
            annotation_bed12,
            treatment_bams, control_bams,
            min_read_overlap, bootstraps,
            log10_transform,
            processes
        )
    write_output_bed(output_bed, results)


if __name__ == '__main__':
    d3pendr()