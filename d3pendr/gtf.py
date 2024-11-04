import re
import numpy as np


def chunk_gtf_records(gtf_records, processes):
    # read the whole gtf file
    gtf_records = list(gtf_records)
    nrecords = len(gtf_records)
    n, r = divmod(nrecords, processes)
    split_points = ([0] + r * [n + 1] + (processes - r) * [n])
    split_points = np.cumsum(split_points)
    for i in range(processes):
        start = split_points[i]
        end = split_points[i + 1]
        yield gtf_records[start: end]


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
            f'from GTF with feature type {gtf_record[2]}'
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
