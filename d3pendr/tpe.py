from collections import defaultdict, Counter
import numpy as np


def relative_tpe(gene_start, gene_end, aln_start, aln_end, strand, read_end):
    # make tpe position relative to gene start rather than 
    # genomic coordinates
    if read_end == '3':
        if strand == '+':
            return aln_end - gene_start
        else:
            return gene_end - aln_start
    elif read_end == '5':
        if strand == '+':
            return aln_start - gene_start
        else:
            return gene_end - aln_end


def intersect(inv_a, inv_b):
    a_start, a_end = inv_a
    b_start, b_end = inv_b
    if a_end < b_start or a_start > b_end:
        return 0
    else:
        s = max(a_start, b_start)
        e = min(a_end, b_end)
        return e - s


def intersect_spliced_invs(invs_a, invs_b):
    score = 0
    invs_a = iter(invs_a)
    invs_b = iter(invs_b)
    a_start, a_end = next(invs_a)
    b_start, b_end = next(invs_b)
    while True:
        if a_end < b_start:
            try:
                a_start, a_end = next(invs_a)
            except StopIteration:
                break
        elif a_start > b_end:
            try:
                b_start, b_end = next(invs_b)
            except StopIteration:
                break
        else:
            score += intersect([a_start, a_end], [b_start, b_end])
            if a_end > b_end:
                try:
                    b_start, b_end = next(invs_b)
                except StopIteration:
                    break
            else:
                try:
                    a_start, a_end = next(invs_a)
                except StopIteration:
                    break
    return score


def get_tpe_distribution(mbam, chrom, gene_start, gene_end, strand,
                         min_read_overlap, read_strand, read_end,
                         paired_end_read, is_bam, as_coverage=False):
    tpe_distribs = []
    nreads = []

    # fetch parsed alignments from bam files, filtering by strand
    if is_bam:

        for sample in mbam.fetch(
                chrom, gene_start, gene_end,
                strand=strand, # the strand that the gene feature is on
                read_strand=read_strand, # the relative orientation of the reads to fetch
                pairs=paired_end_read): # which read end to fetch
            # sample is an iterator of parsed bam file alignments
            sample_distrib = []
            for aln_start, aln_end, aln_strand, aln_invs in sample:
                # calculate the fraction of the read alignment overlapping the
                # annotated gene
                aln_len = sum([e - s for s, e in aln_invs])
                i = intersect_spliced_invs(aln_invs, [(gene_start, gene_end),])

                # only use the read if it overlaps with the reference annotation
                # by at least min_read_overlap fraction of its aligned length
                if i / aln_len >= min_read_overlap:
                    tpe = relative_tpe(
                        gene_start, gene_end,
                        aln_start, aln_end,
                        strand, read_end
                    )
                    sample_distrib.append(tpe)
            if as_coverage:
                ln = gene_end - gene_start
                sample_distrib = np.bincount(sample_distrib, minlength=ln)[:ln]
                n = sum(sample_distrib)
            else:
                n = len(sample_distrib)
            tpe_distribs.append(np.array(sample_distrib))
            nreads.append(n)

    else:
        for sample in mbam.coverage(
                chrom, gene_start, gene_end,
                strand=fetch_strand, normalise=False):
            # sample must be count data
            if np.mod(sample, 1).any():
                raise ValueError('Bigwig values are not discrete')
            sample = sample.astype(np.int)
            # sample is a hist of coverage over the gene from the bigwig file
            # it needs inverting for negative strand genes
            if strand == '-':
                sample = sample[::-1]
            if not as_coverage:
                sample_distrib = np.repeat(np.arange(gene_end - gene_start), sample)
                n = len(sample_distrib)
            else:
                n = sum(sample_distrib)
            tpe_distribs.append(sample_distrib)
            nreads.append(n)

    nreads = np.array(nreads)
    return tpe_distribs, nreads