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
    if strand is None or strand == '.':
        for aln in bam.fetch(*args, **kwargs):
            if pair_filt(aln):
                yield bam_cigar_to_invs(aln)
    elif strand in '+-':
        is_reverse = strand == '-'
        for aln in bam.fetch(*args, **kwargs):
            if is_reverse == aln.is_reverse and pair_filt(aln):
                yield bam_cigar_to_invs(aln)
    else:
        raise ValueError('strand is not one of +-.')


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
            elif strand is None or strand == '.':
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


def flatten_intervals(invs):
    flattened = []
    all_invs = iter(np.sort(invs, axis=0))
    inv_start, inv_end = next(all_invs)
    for start, end in all_invs:
        if start <= inv_end:
            inv_end = max(inv_end, end)
        else:
            flattened.append([inv_start, inv_end])
            inv_start, inv_end = start, end
    if not flattened or flattened[-1] != [inv_start, inv_end]:
        flattened.append([inv_start, inv_end])
    return np.array(flattened)


def filter_terminal_exons(invs, max_intron_size, min_exon_size):
    if len(invs) == 1:
        return invs
    else:
        l_ex = invs[0, 1] - invs[0, 0]
        l_in = invs[1, 0] - invs[0, 1]
        if (l_ex < min_exon_size) or (l_in >= max_intron_size):
            invs = invs[1:]
            if len(invs) == 1:
                return invs
        else:
            r_ex = invs[-1, 1] - invs[-1, 0]
            r_in = invs[-1, 0] - invs[-2, 1]
            if (r_ex < min_exon_size) or (r_in >= max_intron_size):
                invs = invs[:-1]
    return invs


def get_record_range(invs,
                     max_terminal_intron_size=None,
                     min_terminal_exon_size=None,
                     filter_=True):
    invs = flatten_intervals(invs)
    if filter_:
        invs = filter_terminal_exons(
            invs,
            max_terminal_intron_size,
            min_terminal_exon_size,
        )
    return invs[0, 0], invs[-1, 1]


def get_gtf_attribute(gtf_record, attribute):
    try:
        attr = re.search(f'{attribute} "(.+?)";', gtf_record[8]).group(1)
    except AttributeError:
        raise ValueError(
            f'Could not parse attribute {attribute} '
            f'from GTF with feature type {record[2]}'
        )
    return attr


def gtf_iterator(gtf_fn,
                 extend_gene_five_prime=0,
                 use_5utr=False,
                 extend_gene_three_prime=0,
                 by_locus=True,
                 max_terminal_intron_size=100_000,
                 min_terminal_exon_size=20):
    gtf_records = {}
    if by_locus:
        gene_to_locus_mapping = {}
    with open(gtf_fn) as gtf:
        for i, record in enumerate(gtf):
            record = record.split('\t')
            chrom, _, feat_type, start, end, _, strand = record[:7]
            start = int(start) - 1
            end = int(end)
            if feat_type == 'transcript' and by_locus:
                locus_id = get_gtf_attribute(record, 'locus')
                gene_id = get_gtf_attribute(record, 'gene_id')
                gene_to_locus_mapping[gene_id] = locus_id
            elif feat_type == 'CDS' or feat_type == 'exon':
                gene_id = get_gtf_attribute(record, 'gene_id')
                idx = (chrom, gene_id, strand)
                if idx not in gtf_records:
                    gtf_records[idx] = {}
                if feat_type not in gtf_records[idx]:
                    gtf_records[idx][feat_type] = []
                gtf_records[idx][feat_type].append((start, end))

    if by_locus:
        # regroup gene invs by locus id:
        gtf_records_by_locus = {}
        for (chrom, gene_id, strand), feat_invs in gtf_records.items():
            locus_id = gene_to_locus_mapping[gene_id]
            new_idx = (chrom, locus_id, strand)
            if new_idx not in gtf_records_by_locus:
                gtf_records_by_locus[new_idx] = {}
            for feat_type, invs in feat_invs.items():
                if feat_type not in gtf_records_by_locus[new_idx]:
                    gtf_records_by_locus[new_idx][feat_type] = []
                gtf_records_by_locus[new_idx][feat_type] += invs
        gtf_records = gtf_records_by_locus

    # once whole file is parsed yield the intervals
    for (chrom, gene_id, strand), feat_invs in gtf_records.items():
        exon_start, exon_end = get_record_range(
            feat_invs['exon'],
            max_terminal_intron_size,
            min_terminal_exon_size,
        )
        try:
            cds_start, cds_end = get_record_range(
                feat_invs['CDS'],
                filter_=False,
            )
        except KeyError:
            # non-coding RNA
            cds_start, cds_end = exon_start, exon_end

        # remove region corresponding to 5'UTR if necessary
        if use_5utr:
            gene_start = exon_start
            gene_end = exon_end
        else:
            gene_start = cds_start if strand == '+' else exon_start
            gene_end = exon_end if strand == '+' else cds_end
        

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