import dataclasses
import click

from .output import write_wass_test_output_bed
from .d3pendr import run_d3pendr


@dataclasses.dataclass
class BAMOpts:
    treatment_fns: str
    control_fns: str
    read_strand: str
    read_end: str
    paired_end_read: str
    min_read_overlap: float
    min_reads_per_rep: int

@dataclasses.dataclass
class GTFOpts:
    annotation_gtf_fn: str
    extend_gene_five_prime: int
    extend_gene_three_prime: int
    use_locus_tag: bool
    allow_overlap_next_gene: bool


@dataclasses.dataclass
class StatOpts:
    bootstraps: int
    threshold: float
    wass_fit_gamma: bool
    silhouette_test: bool
    sil_fit_skewnorm: bool

        
@dataclasses.dataclass(init=False)
class D3pendrOpts:
    bam: BAMOpts
    gtf: GTFOpts
    stats: StatOpts
    output_bed: str
    processes: int
    random_seed: int = None


    def __init__(self, **kwargs):
        _nested_fields = {}
        _toplevel_fields = set()
        for field in dataclasses.fields(self):
            if dataclasses.is_dataclass(field.type):
                _nested_fields[(field.name, field.type)] = set()
                for subfield in dataclasses.fields(field.type):
                    _nested_fields[(field.name, field.type)].add(subfield.name)
            else:
                _toplevel_fields.add(field.name)
        # assign values to sub-dataclasses
        for (name, subclass), fields in _nested_fields.items():
            subclass_kwargs = {}
            for f in fields:
                subclass_kwargs[f] = kwargs[f]
            setattr(self, name, subclass(**subclass_kwargs))
        # assign toplevel fields
        for name in _toplevel_fields:
            setattr(self, name, kwargs[name])


def make_dataclass_decorator(dc):
    def _dataclass_decorator(cmd):
        @click.pass_context
        def _make_dataclass(ctx, **kwargs):
            return ctx.invoke(cmd, dc(**kwargs))
        return _make_dataclass
    return _dataclass_decorator


@click.command()
@click.option('-t', '--treatment-fns', required=True, multiple=True)
@click.option('-c', '--control-fns', required=True, multiple=True)
@click.option('-o', '--output-bed', required=True)
@click.option('-a', '--annotation-gtf-fn', required=True)
@click.option('--read-strand', type=click.Choice(['same', 'opposite', 'unstranded']), default='same')
@click.option('--read-end', type=click.Choice(['3', '5']), default='3')
@click.option('--paired-end-read', type=click.Choice(['1', '2', 'both', 'single']), default='single')
@click.option('--min-read-overlap', default=0.2)
@click.option('--min-reads-per-rep', default=5)
@click.option('--extend-gene-five-prime', default=0)
@click.option('--extend-gene-three-prime', default=0)
@click.option('--use-locus-tag/--use-gene-id-tag', default=False)
@click.option('--allow-overlap-next-gene/--no-gene-overlap', default=False)
@click.option('--bootstraps', default=999)
@click.option('--threshold', default=0.05)
@click.option('--wass-fit-gamma/--no-fit-gamma', default=True)
@click.option('--silhouette-test/--no-silhouette-test', default=True)
@click.option('--sil-fit-skewnorm/--no-fit-skewnorm', default=True)
@click.option('-p', '--processes', default=4)
@click.option('-r', '--random-seed', default=None)
@make_dataclass_decorator(D3pendrOpts)
def d3pendr(opts):
    '''
    d3pendr: Differential 3' End analysis of Nanopore Direct RNAs

    Identifies differential polyadenylation events in either an
    annotation dependent manner, using
    a permutation test of Wasserstein distance between 3' end 
    distributions and a KS test.

    Outputs bed6 format with extra columns.
    '''

    if opts.random_seed is None:
        opts.random_seed = abs(hash(opts.gtf.annotation_gtf_fn))

    results = run_d3pendr(opts)
    write_wass_test_output_bed(opts.output_bed, results)


if __name__ == '__main__':
    d3pendr()