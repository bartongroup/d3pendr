def write_wass_test_output_bed(output_bed_fn, results):
    with open(output_bed_fn, 'w') as bed:
        for (
            chrom, start, end, gene_id, score, strand,
            wass_dist, wass_dir, wass_pval, wass_fdr,
            sil_score, sil_pval, sil_fdr,
            nreads_cntrl, nreads_treat
        ) in results.itertuples(index=False):
            record = (f'{chrom}\t{int(start):d}\t{int(end):d}\t'
                      f'{gene_id}\t{int(score):d}\t{strand}\t'
                      f'{wass_dist:.1f}\t{wass_dir:.1f}\t'
                      f'{wass_pval:.3g}\t{wass_fdr:.3g}\t'
                      f'{sil_score:.3f}\t{sil_pval:.3g}\t{sil_fdr:.3g}\t'
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
