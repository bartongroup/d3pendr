import re
import numpy as np
import dataclasses


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


def get_gtf_attribute(gtf_record, attribute):
    try:
        attr = re.search(f'{attribute} "(.+?)";', gtf_record[8]).group(1)
    except AttributeError:
        raise ValueError(
            f'Could not parse attribute {attribute} '
            f'from GTF with feature type {record[2]}'
        )
    return attr


@dataclasses.dataclass
class GTFRecord:
    '''Class for handling results of wasserstein_test and wasserstein_silhouette_test'''
    chrom: str
    start: int
    end: int
    locus_id: str
    strand: str

    @property
    def is_reverse(self):
        return self.strand == '-'


import re
import numpy as np
import dataclasses


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


def get_gtf_attribute(gtf_record, attribute):
    try:
        attr = re.search(f'{attribute} "(.+?)";', gtf_record[8]).group(1)
    except AttributeError:
        raise ValueError(
            f'Could not parse attribute {attribute} '
            f'from GTF with feature type {record[2]}'
        )
    return attr


@dataclasses.dataclass
class GTFRecord:
    '''Class for handling results of wasserstein_test and wasserstein_silhouette_test'''
    chrom: str
    start: int
    end: int
    locus_id: str
    strand: str

    @property
    def is_reverse(self):
        return self.strand == '-'


class GTFGeneBoundaryIterator:

    def __init__(self, annotation_gtf_fn,
                 extend_gene_five_prime=0,
                 extend_gene_three_prime=0,
                 use_locus_tag=False,
                 allow_overlap_next_gene=False):
        self._fn = annotation_gtf_fn
        self._use_locus_tag = use_locus_tag
        self._extend_5p = extend_gene_five_prime
        self._extend_3p = extend_gene_three_prime
        self._allow_overlap = allow_overlap_next_gene
        self._gtf_records = self._parse()

    def _parse(self):
        gtf_records = {}
        if self._use_locus_tag:
            gene_to_locus_mapping = {}
        with open(self._fn) as gtf:
            for record in gtf:
                record = record.split('\t')
                chrom, _, feat_type, start, end, _, strand = record[:7]
                start = int(start) - 1
                end = int(end)
                if feat_type == 'transcript' and use_locus_tag:
                    locus_id = get_gtf_attribute(record, 'locus')
                    gene_id = get_gtf_attribute(record, 'gene_id')
                    gene_to_locus_mapping[gene_id] = locus_id
                elif feat_type == 'exon':
                    gene_id = get_gtf_attribute(record, 'gene_id')
                    if gene_id not in gtf_records:
                        gtf_records[gene_id] = GTFRecord(
                            chrom, start, end, gene_id, strand
                        )
                    else:
                        gtf_records[gene_id].start = min(
                            start, gtf_records[gene_id].start
                        )
                        gtf_records[gene_id].end = max(
                            end, gtf_records[gene_id].end
                        )

        if self._use_locus_tag:
            gtf_records_by_locus = {}
            for gene_id, record in gtf_records.items():
                locus_id = gene_to_locus_mapping[gene_id]
                if gene_id not in gtf_records_by_locus:
                    gtf_records_by_locus[locus_id] = record
                else:
                    gtf_records_by_locus[locus_id].start = min(
                        record.start, gtf_records_by_locus[locus_id].start
                    )
                    gtf_records_by_locus[locus_id].end = max(
                        record.end, gtf_records_by_locus[locus_id].end
                    )
            gtf_records = gtf_records_by_locus
        
        # reorganise by chrom and strand
        gtf_records_by_chom = {}
        for record in gtf_records.values():
            idx = (record.chrom, record.strand)
            if idx not in gtf_records_by_chom:
                gtf_records_by_chom[idx] = []
            gtf_records_by_chom[idx].append(record)

        for records in gtf_records_by_chom.values():
            records.sort(key=lambda r: (r.start, r.end))
        self._records = gtf_records_by_chom

    def __iter__(self):
        for records in self._records.values():
            for i, r in enumerate(records):
                # add extensions to 3' and 5' ends
                if r.strand == '+':
                    extend_start, extend_end = self._extend_5p, self._extend_3p
                else:
                    extend_start, extend_end = self._extend_3p, self._extend_5p

                if not self._allow_overlap:
                    if i != 0:
                        extend_start = min(
                            extend_start,
                            r.start - records[i - 1].end
                        )

                    try:
                        extend_end = min(
                            extend_end,
                            records[i + 1].start - r.end
                        )
                    except IndexError:
                        pass

                r.start = max(0, r.start - extend_start)
                r.end = r.end + extend_end
                yield r

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass