import re
from functools import lru_cache

import pysam
import click


def pair_filter(filt_type):
    if filt_type == '1':
        def _pair_filter(aln):
            return aln.is_read1
    elif filt_type == '2':
        def _pair_filter(aln):
            return aln.is_read2
    return _pair_filter


def get_polya_test_coords(aln, end, range_):
    chrom = aln.reference_name
    leftmap = aln.reference_start
    rightmap = aln.reference_end
    strand = '+-'[aln.is_reverse]
    if strand == '+':
        if end == '5':
            start = max(leftmap - range_[0], 0)
            end = leftmap + range_[1]
        else:
            start = rightmap - range_[1]
            end = rightmap + range_[0]
    elif strand == '-':
        if end == '5':
            start = rightmap - range_[1]
            end = rightmap + range_[0]
        else:
            start = max(leftmap - range_[0], 0)
            end = leftmap + range_[1]
    return chrom, start, end, strand


def filter_internal_priming(aln_stream,
                            fasta,
                            is_paired=False,
                            polyt_read='1',
                            search_end='5',
                            search_range=(10, 5),
                            min_as_discontinuous=7,
                            min_a_homopolymer=6,
                            min_purines=8,
                            min_as_in_purine_rich=4,
                            min_purine_homopolymer=7,
                            min_as_in_purine_homopolymer=4):
    search_len = sum(search_range)
    if is_paired:
        pair_filt = pair_filter(polyt_read)
    else:
        pair_filt = lambda aln: True

    RC = str.maketrans('ACGTRYSWN', 'TGCAYRSWN')

    @lru_cache(maxsize=128)
    def fasta_fetch(chrom, start, end, strand):
        seq = fasta.fetch(chrom, start, end).upper()
        if strand == '+':
            seq = seq.translate(RC)[::-1]
        return seq

    for aln in aln_stream:
        if pair_filt(aln):
            chrom, start, end, strand = get_polya_test_coords(
                aln, search_end, search_range
            )
            search_seq = fasta_fetch(chrom, start, end, strand)

            a_count = search_seq.count('A')
            g_count = search_seq.count('G')
            purine_count = a_count + g_count

            if a_count >= min_as_discontinuous:
                continue
            elif (purine_count >= min_purines) and \
                 (a_count >= min_as_in_purine_rich):
                continue
            elif re.search(f'A{{{min_a_homopolymer},}}', search_seq):
                continue
            else:
                # final test, if using python 3.8 could use assignment expression
                purine_homopolymers = re.findall(
                    f'[AG]{{{min_as_in_purine_homopolymer},}}',
                    search_seq
                )
                for hp in purine_homopolymers:
                    if hp.count('A') >= min_as_in_purine_homopolymer:
                        break
                else:
                    yield aln


@click.command()
@click.option('-b', '--bam-fn', required=True)
@click.option('-o', '--output-bam-fn', required=True)
@click.option('-f', '--fasta-fn', required=True)
@click.option('--paired/--single', default=False)
@click.option('--polyt-primed-read-in-pair', type=click.Choice(['1', '2']), default='1')
@click.option('--polyt-primed-read-end', type=click.Choice(['5', '3']), default='5')
@click.option('--search-range', nargs=2, default=(10, 5))
@click.option('--min-a-discontinuous', default=7)
@click.option('--min-a-homopolymer', default=6)
@click.option('--min-purine-discontinuous', default=9)
@click.option('--min-a-in-purine-rich', default=4)
@click.option('--min-purine-homopolymer', default=7)
@click.option('--min-a-in-purine-homopolymer', default=4)
def main(bam_fn, output_bam_fn, fasta_fn, paired,
         polyt_primed_read_in_pair, polyt_primed_read_end, search_range,
         min_a_discontinuous, min_a_homopolymer,
         min_purine_discontinuous, min_a_in_purine_rich,
         min_purine_homopolymer, min_a_in_purine_homopolymer):
    with pysam.AlignmentFile(bam_fn) as inbam, pysam.FastaFile(fasta_fn) as fasta:

        with pysam.AlignmentFile(output_bam_fn, 'wb', template=inbam) as outbam:

            filtered_alns = filter_internal_priming(
                inbam.fetch(),
                fasta,
                paired,
                polyt_primed_read_in_pair,
                polyt_primed_read_end,
                search_range,
                min_a_discontinuous,
                min_a_homopolymer,
                min_purine_discontinuous,
                min_a_in_purine_rich,
                min_purine_homopolymer,
                min_a_in_purine_homopolymer,
            )
            for aln in filtered_alns:
                outbam.write(aln)


if __name__ == '__main__':
    main()