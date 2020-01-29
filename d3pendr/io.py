import heapq
from operator import itemgetter
import numpy as np
import pysam


def bam_cigar_to_invs(aln):
    invs = []
    start = aln.reference_start
    end = aln.reference_end
    strand = '-' if aln.is_reverse else '+'
    left = start
    right = left
    has_ins = False
    for op, ln in aln.cigar:
        if op in (1, 4, 5):
            # does not consume reference
            continue
        elif op in (0, 2, 7, 8):
            # consume reference but do not add to invs yet
            right += ln
        elif op == 3:
            invs.append([left, right])
            left = right + ln
            right = left
    if right > left:
        invs.append([left, right])
    assert invs[0][0] == start
    assert invs[-1][1] == end
    return start, end, strand, np.array(invs)


def bam_query_iterator(bam, *args, **kwargs):
    strand = kwargs.pop('strand', None)
    if strand is None:
        for aln in bam.fetch(*args, **kwargs):
            yield bam_cigar_to_invs(aln)
    elif strand in '+-':
        is_reverse = strand == '-'
        for aln in bam.fetch(*args, **kwargs):
            if is_reverse == aln.is_reverse:
                yield bam_cigar_to_invs(aln)
    else:
        raise ValueError('strand is not one of +-')
                


class MultiBamParser(object):

    def __init__(self, bam_fns, *, sort_fetch=False):
        self.bam_handles = {
            bam_fn: pysam.AlignmentFile(bam_fn) for bam_fn in bam_fns
        }
        self.closed = False
        self._sort_fetch = sort_fetch

    def fetch(self, *args, **kwargs):
        queries = [
            bam_query_iterator(bam, *args, **kwargs)
            for bam in self.bam_handles.values()
        ]
        if self._sort_fetch:
            yield from heapq.merge(*queries, key=itemgetter(1))
        else:
            for q in queries:
                yield from q

    def close(self):
        for bam in self.bam_handles.values():
            bam.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def parse_exons(record):
    start = int(record[1])
    end = int(record[2])
    exstarts = np.fromstring(record[11], sep=',') + start
    exends = exstarts + np.fromstring(record[10], sep=',')
    exons = np.dstack([exstarts, exends])[0]
    return start, end, exons


def parse_bed_record(record):
    chrom = record[0]
    strand = record[5]
    gene_id = record[3]
    start, end, invs = parse_exons(record)
    return chrom, start, end, gene_id, strand, invs


def bed12_iterator(bed_fn):
    with open(bed_fn) as bed:
        for record in bed:
            record = record.split()
            yield parse_bed_record(record)


def write_output_bed(output_bed_fn, results):
    with open(output_bed_fn, 'w') as bed:
        for (
            chrom, start, end, gene_id, score, strand,
            wass_dist, wass_dir, wass_pval, wass_fdr,
            ks_stat, ks_pval, ks_fdr, hm_pval, hm_fdr,
            nreads_cntrl, nreads_treat
        ) in results.itertuples(index=False):

            record = (f'{chrom}\t{int(start):d}\t{int(end):d}\t'
                      f'{gene_id}\t{int(score):d}\t{strand}\t'
                      f'{wass_dist:.1f}\t{wass_dir:.1f}\t'
                      f'{wass_pval:.3g}\t{wass_fdr:.3g}\t'
                      f'{ks_stat:.2f}\t{ks_pval:.3g}\t{ks_fdr:.3g}\t'
                      f'{hm_pval:.3g}\t{hm_fdr:.3g}\t'
                      f'{nreads_cntrl:d}\t{nreads_treat:d}\n')
            bed.write(record)