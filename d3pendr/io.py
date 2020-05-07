import numpy as np
import pysam
import pyBigWig as pybw


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


def pair_filter(filt_type):
    if filt_type == 'both':
        def _pair_filter(aln):
            return True
    elif filt_type == '1':
        def _pair_filter(aln):
            return aln.is_read1
    elif filt_type == '2':
        def _pair_filter(aln):
            return aln.is_read2
    return _pair_filter


def bam_query_iterator(bam, *args, **kwargs):
    strand = kwargs.pop('strand', None)
    pairs = kwargs.pop('pairs', 'both')
    pair_filt = pair_filter(pairs)
    if strand is None:
        for aln in bam.fetch(*args, **kwargs):
            if pair_filt(aln):
                yield bam_cigar_to_invs(aln)
    elif strand in '+-':
        is_reverse = strand == '-'
        for aln in bam.fetch(*args, **kwargs):
            if is_reverse == aln.is_reverse and pair_filt(aln):
                yield bam_cigar_to_invs(aln)
    else:
        raise ValueError('strand is not one of +-')


class MultiParser(object):

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class MultiBamParser(MultiParser):

    def __init__(self, bam_fns):
        self.handles = {
            bam_fn: pysam.AlignmentFile(bam_fn) for bam_fn in bam_fns
        }
        self.closed = False

    def fetch(self, *args, **kwargs):
        queries = [
            bam_query_iterator(bam, *args, **kwargs)
            for bam in self.handles.values()
        ]
        return queries

    def close(self):
        for bam in self.handles.values():
            bam.close()


class MultiBigWigParser(MultiParser):

    def __init__(self, bw_fns):
        # attempt to infer if stranded, each bw_fn should be comma separated list
        bw_fns = [tuple(fn.split(',')) for fn in bw_fns]
        if all([len(fn) == 2 for fn in bw_fns]):
            stranded = True
        elif all([len(fn) == 1 for fn in bw_fns]):
            stranded = False
        else:
            raise ValueError('Please provide either single bw files or comma separated pos,neg bw files')
        if stranded:
            self.handles = {
                bw_fn: (pybw.open(bw_fn[0]), pybw.open(bw_fn[1])) for bw_fn in bw_fns
            }
        else:
            self.handles = {
                bw_fn: (pybw.open(bw_fn[0]),) for bw_fn in bw_fns
            }
        self.closed = False
        self.stranded = stranded

    def fetch(self, chrom, start, end, strand=None):
        if strand is not None and not self.stranded:
            raise ValueError('cannot specify strand on unstranded bigwigs')
        if self.stranded:
            if strand == '+':
                queries = [
                    pos_bw.values(chrom, start, end, numpy=True)
                    for pos_bw, neg_bw in self.handles.values()
                ]
            elif strand == '-':
                queries = [
                    neg_bw.values(chrom, start, end, numpy=True)
                    for pos_bw, neg_bw in self.handles.values()
                ]
            elif strand is None:
                queries = [
                    np.nansum([
                        pos_bw.values(chrom, start, end, numpy=True),
                        neg_bw.values(chrom, start, end, numpy=True)
                    ], axis=0)
                    for pos_bw, neg_bw in self.handles.values()
                ]
        else:
            queries = [
                bw[0].values(chrom, start, end, numpy=True)
                for bw in self.handles
            ]
        return queries

    def close(self):
        for bws in self.handles.values():
            for bw in bws:
                bw.close()


def bam_or_bw(fn):
    if fn.endswith('.bam') or fn.endswith('.sam'):
        return True
    elif fn.endswith('.bw') or fn.lower().endswith('.bigwig'):
        return False
    else:
        raise ValueError('files must be bam, sam, or bigwig format')


def parse_inv(record, use_5utr):
    chrom = record[0]
    start = int(record[1])
    end = int(record[2])
    gene_id = record[3]
    strand = record[5]
    if not use_5utr:
        cds_start = int(record[6])
        cds_end = int(record[7])
        if cds_start != cds_end:
            # not a protein coding gene
            if strand == '+':
                start = cds_start
            else:
                end = cds_end
    return chrom, start, end, strand, gene_id


def parse_bed_record(record, extend_gene_five_prime=0, use_5utr=True, extend_gene_three_prime=0):
    chrom, start, end, strand, gene_id = parse_inv(record, use_5utr)
    if extend_gene_five_prime:
        if strand == '+':
            start = max(0, start - extend_gene_five_prime)
        else:
            end += extend_gene_five_prime
    if extend_gene_three_prime:
        if strand == '+':
            end += extend_gene_three_prime
        else:
            start = max(0, start - extend_gene_three_prime)
    return chrom, start, end, gene_id, strand


def bed12_iterator(bed_fn, extend_gene_five_prime, ignore_5utr, extend_gene_three_prime):
    with open(bed_fn) as bed:
        for record in bed:
            record = record.split()
            yield parse_bed_record(
                record, extend_gene_five_prime, ignore_5utr, extend_gene_three_prime
            )


def write_output_bed(output_bed_fn, results):
    with open(output_bed_fn, 'w') as bed:
        for (
            chrom, start, end, gene_id, score, strand,
            wass_dist, wass_dir, wass_pval, wass_fdr,
            nreads_cntrl, nreads_treat
        ) in results.itertuples(index=False):
            record = (f'{chrom}\t{int(start):d}\t{int(end):d}\t'
                      f'{gene_id}\t{int(score):d}\t{strand}\t'
                      f'{wass_dist:.1f}\t{wass_dir:.1f}\t'
                      f'{wass_pval:.3g}\t{wass_fdr:.3g}\t'
                      f'{sum(nreads_cntrl):d}\t{sum(nreads_treat):d}\n')
            bed.write(record)


def write_apa_site_bed(output_bed_fn, results):
    with open(output_bed_fn, 'w') as bed:
        for (
            chrom, start, end, gene_id, strand,
            cntrl_count, treat_count,
            cntrl_frac, treat_frac,
            relative_change
        ) in results.itertuples(index=False):
            direction = int(relative_change > 0)
            record = (f'{chrom}\t{start:d}\t{end:d}\t'
                      f'{gene_id}\t{direction:d}\t{strand}\t'
                      f'{cntrl_count:d}\t{treat_count:d}\t'
                      f'{cntrl_frac:.2f}\t{treat_frac:.2f}\t'
                      f'{relative_change:.2f}\n')
            bed.write(record)