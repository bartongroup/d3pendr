import numpy as np

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
    if filt_type == 'both' or filt_type == 'single':
        def _pair_filter(aln):
            return True
    elif filt_type == '1':
        def _pair_filter(aln):
            return aln.is_read1
    elif filt_type == '2':
        def _pair_filter(aln):
            return aln.is_read2
    return _pair_filter


def strand_filter(filt_type, pair_type):
    if pair_type != "both":
        # if we are just using either read1 or 2,
        # assume that "same" or "opposite" refer
        # to the orientation of that read to the
        # strand
        if filt_type == 'same':
            def _strand_filter(aln, strand):
                is_reverse = strand == '-'
                return aln.is_reverse == is_reverse
        if filt_type == 'opposite':
            def _strand_filter(aln, strand):
                is_reverse = strand == '-'
                return aln.is_reverse != is_reverse
    else:
        # in paired ended data, assume that reads1&2 are
        # on the opposite strand from each other and
        # "same" or "opposite" is referring to read1
        # i.e. in paired mode, "same" returns read1 from the same
        # strand and read2 from the opposite strand
        if filt_type == 'same':
            def _strand_filter(aln, strand):
                is_reverse = strand == '-'
                strand_same = aln.is_reverse == is_reverse
                return aln.is_read1 == strand_same
        if filt_type == 'opposite':
            def _strand_filter(aln, strand):
                is_reverse = strand == '-'
                strand_same = aln.is_reverse == is_reverse
                return aln.is_read1 != strand_same
    return _strand_filter



def bam_query_iterator(bam, *args, **kwargs):
    strand = kwargs.pop('strand', None) # + or -
    orient = kwargs.pop('read_strand', 'same') # same or opposite
    pairs = kwargs.pop('pairs', 'single') # both, single, "1" or "2"
    strand_filt = strand_filter(orient, pairs)
    pair_filt = pair_filter(pairs)
    if strand is None or strand == '.':
        for aln in bam.fetch(*args, **kwargs):
            if pair_filt(aln):
                yield bam_cigar_to_invs(aln)
    elif strand in '+-':
        for aln in bam.fetch(*args, **kwargs):
            if strand_filt(aln, strand) and pair_filt(aln):
                yield bam_cigar_to_invs(aln)
    else:
        raise ValueError('strand is not one of +-.')


def add_read_to_cov(cov, aln, offset):
    *_, invs = bam_cigar_to_invs(aln)
    for s, e in invs:
        s = max(s - offset, 0)
        e = e - offset
        cov[s: e] += 1
    return cov


def per_base_bam_coverage(bam, chrom, start, end, **kwargs):
    strand = kwargs.pop('strand', None)
    orient = kwargs.pop('read_strand', 'same')
    pairs = kwargs.pop('pairs', 'both')
    strand_filt = strand_filter(orient, pairs)
    pair_filt = pair_filter(pairs)
    cov = np.zeros(end - start, dtype='int64')
    if strand is None or strand == '.':
        for aln in bam.fetch(chrom, start, end, **kwargs):
            if pair_filt(aln):
                cov = add_read_to_cov(cov, aln, start)
    elif strand in '+-':
        is_reverse = strand == '-'
        for aln in bam.fetch(chrom, start, end, **kwargs):
            if strand_filt(aln, strand) and pair_filt(aln):
                cov = add_read_to_cov(cov, aln, start)
    else:
        raise ValueError('strand is not one of +-.')
    return cov