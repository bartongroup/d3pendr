import re
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


def gtf_iterator(gtf_fn, extend_gene_five_prime, ignore_5utr, extend_gene_three_prime):
    gtf_records = {}
    with open(gtf_fn) as gtf:
        for i, record in enumerate(gtf):
            record = record.split('\t')
            feat_type = record[2]
            if feat_type == 'CDS' or feat_type == 'exon':
                try:
                    gene_id = re.search('gene_id "(.+?)";', record[8]).group(1)
                except AttributeError:
                    raise ValueError(f'Could not parse gene_id from GTF line {i}')
                if gene_id not in gtf_records:
                    gtf_records[gene_id] = {
                        'chrom': record[0],
                        'strand': record[6]
                    }
                start = int(record[3]) - 1
                end = int(record[4])
                if feat_type not in gtf_records[gene_id]:
                    gtf_records[gene_id][feat_type] = (start, end)
                else:
                    curr_range = gtf_records[gene_id][feat_type]
                    new_range = (
                        min(curr_range[0], start),
                        max(curr_range[1], end),
                    )
                    gtf_records[gene_id][feat_type] = new_range

    # once whole file is parsed yield the intervals
    for gene_id, gene_info in gtf_records.items():
        chrom = gene_info['chrom']
        strand = gene_info['strand']
        exon_start, exon_end = gene_info['exon']
        try:
            cds_start, cds_end = gene_info['CDS']
        except KeyError:
            # non-coding RNA
            cds_start, cds_end = gene_info['exon']

        # remove region corresponding to 5'UTR if necessary
        if ignore_5utr:
            gene_start = cds_start if strand == '+' else exon_start
            gene_end = exon_end if strand == '+' else cds_end
        else:
            gene_start = exon_start
            gene_end = exon_end

        # add extensions to 3' and 5' ends
        start_ext, end_ext = extend_gene_five_prime, extend_gene_three_prime
        if strand == '-':
            start_ext, end_ext = end_ext, start_ext
        gene_start = max(0, gene_start - start_ext)
        gene_end = gene_end + end_ext

        yield chrom, gene_start, gene_end, gene_id, strand


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