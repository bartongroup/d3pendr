import numpy as np

import pysam
import pyBigWig as pybw

from .bam import bam_query_iterator, per_base_bam_coverage

def _bam_or_bw(fn):
    if fn.endswith('.bam') or fn.endswith('.sam'):
        return True
    elif fn.endswith('.bw') or fn.lower().endswith('.bigwig'):
        return False
    else:
        raise ValueError('files must be bam, sam, or bigwig format')


def bam_or_bw(*fns):
    filetypes = [_bam_or_bw(fn) for fn in fns]
    assert all([t == filetypes[0] for t in filetypes])
    filetype = filetypes[0]
    return filetype


class MultiParser(object):

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class MultiBamParser(MultiParser):

    def __init__(self, bam_fns):
        self.handles = [
            pysam.AlignmentFile(bam_fn) for bam_fn in bam_fns
        ]
        self._calc_norm_factors(self.handles)
        self.closed = False

    def fetch(self, *args, **kwargs):
        queries = [
            bam_query_iterator(bam, *args, **kwargs)
            for bam in self.handles
        ]
        return queries

    def coverage(self, chrom, start, end, **kwargs):
        normalise = kwargs.pop('normalise', False)
        coverage = [
            per_base_bam_coverage(bam, chrom, start, end, **kwargs)
            for bam in self.handles
        ]
        coverage = np.array(coverage)

        if normalise:
            coverage *= self.norm_factors
        
        return coverage

    def _calc_norm_factors(self, bams):
        p = np.reshape([b.mapped for b in bams], newshape=(-1, 1))
        n = p.mean()
        self.norm_factors = n / p

    def close(self):
        for bam in self.handles:
            bam.close()


class MultiBigWigParser(MultiParser):

    def __init__(self, bw_fns):
        # attempt to infer if stranded, each bw_fn should be 
        # comma separated list
        bw_fns = [tuple(fn.split(',')) for fn in bw_fns]
        if all([len(fn) == 2 for fn in bw_fns]):
            stranded = True
        elif all([len(fn) == 1 for fn in bw_fns]):
            stranded = False
        else:
            raise ValueError('Please provide either single bw files '
                             'or comma separated pos,neg bw files')
        if stranded:
            self.handles = [
                (pybw.open(bw_fn[0]), pybw.open(bw_fn[1]))
                for bw_fn in bw_fns
            ]
        else:
            self.handles = [
                (pybw.open(bw_fn[0]),) for bw_fn in bw_fns
            ]
        self._calc_norm_factors(self.handles)
        self.closed = False
        self.stranded = stranded

    def fetch(self, *args, **kwargs):
        raise NotImplemented('fetch not possible for bigwig')

    def coverage(self, chrom, start, end, strand=None, normalise=False, **kwargs):
        if strand is not None and not self.stranded:
            raise ValueError('cannot specify strand on unstranded bigwigs')
        if self.stranded:
            if strand == '+':
                coverage = [
                    pos_bw.values(chrom, start, end, numpy=True)
                    for pos_bw, neg_bw in self.handles
                ]
            elif strand == '-':
                coverage = [
                    neg_bw.values(chrom, start, end, numpy=True)
                    for pos_bw, neg_bw in self.handles
                ]
            elif strand is None or strand == '.':
                coverage = [
                    np.nansum([
                        pos_bw.values(chrom, start, end, numpy=True),
                        neg_bw.values(chrom, start, end, numpy=True)
                    ], axis=0)
                    for pos_bw, neg_bw in self.handles
                ]
        else:
            coverage = [
                bw[0].values(chrom, start, end, numpy=True)
                for bw in self.handles
            ]
        coverage = np.array(coverage)
        coverage[np.isnan(coverage)] = 0

        if normalise:
            coverage *= self.norm_factors
        
        return coverage

    def _calc_norm_factors(self, bws):
        p = np.reshape(
            [sum([s.header()['sumData'] for s in b]) for b in bws],
            newshape=(-1, 1)
        )
        n = p.mean()
        self.norm_factors = n / p

    def close(self):
        for bws in self.handles:
            for bw in bws:
                bw.close()