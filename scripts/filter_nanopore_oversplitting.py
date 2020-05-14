import numpy as np
import pandas as pd
import pysam
import click


def get_read_mapping_locs(bam_fn):
    mapping = []
    with pysam.AlignmentFile(bam_fn) as bam:
        for aln in bam.fetch():
            mapping.append([
                aln.query_name,
                aln.reference_name,
                aln.reference_start,
                aln.reference_end,
                ['+', '-'][aln.is_reverse],
            ])
    return pd.DataFrame(
        mapping,
        columns=['read_id', 'chrom', 'genomic_start', 'genomic_end', 'strand']
    )


def get_dist_to_prev(read):
    if not read.strand_same_as_prev or not read.chrom_same_as_prev:
        return np.nan
    if read.strand == '+':
        return read.genomic_start_prev_read - read.genomic_end
    elif read.strand == '-':
        return read.genomic_start - read.genomic_end_prev_read


def estimate_oversplitting(ss_fn, bam_fn):
    ss = pd.read_csv(ss_fn, sep='\t', usecols=['read_id', 'channel', 'start_time'])
    ss = ss.sort_values(['channel', 'start_time'])
    ss['prev_read_id'] = ss.groupby('channel').read_id.shift(1)
    mapped_loc = get_read_mapping_locs(bam_fn)
    ss = ss.merge(mapped_loc, on='read_id', how='left')
    ss = ss.merge(
        mapped_loc,
        left_on='prev_read_id',
        right_on='read_id',
        suffixes=('', '_prev_read'),
        how='left'
    ).drop('read_id_prev_read', axis=1)
    ss['chrom_same_as_prev'] = ss.chrom == ss.chrom_prev_read
    ss['strand_same_as_prev'] = ss.strand == ss.strand_prev_read
    ss['genomic_dist_to_prev'] = ss.apply(get_dist_to_prev, axis=1)
    return ss


def filter_bam(input_bam_fn, output_bam_fn, blacklist_read_ids):
    with pysam.AlignmentFile(input_bam_fn) as bam:
        with pysam.AlignmentFile(output_bam_fn, 'wb', template=bam) as output_bam:
            for aln in bam.fetch():
                if not aln.query_name in blacklist_read_ids:
                    output_bam.write(aln)


@click.command()
@click.option('-b', '--bam-fn', required=True)
@click.option('-s', '--sequencing-summary-fn', required=True)
@click.option('-o', '--output-bam-fn', required=True)
@click.option('-d', '--max-distance', default=1000)
@click.option('-d', '--min-distance', default=-10)
def cli(bam_fn, sequencing_summary_fn, output_bam_fn,
        max_distance, min_distance):
    oversplit_reads = estimate_oversplitting(sequencing_summary_fn, bam_fn)
    oversplit_reads = oversplit_reads.dropna(subset=['genomic_dist_to_prev'])
    oversplit_reads = oversplit_reads.query(
        f'{min_distance} <= genomic_dist_to_prev <= {max_distance}'
    )
    # for the purpose of estimating changes in 3' ends, we want to filter out
    # the upstream read as this is the one with the erroneous 3' end
    # upstream read is the one sequenced later (reads are sequenced 3' to 5')
    oversplit_read_ids = set(oversplit_reads.read_id)
    filter_bam(bam_fn, output_bam_fn, oversplit_read_ids)
    pysam.index(output_bam_fn)


if __name__ == '__main__':
    cli()