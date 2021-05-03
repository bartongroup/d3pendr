import numpy as np

def format_bed_record(row):
    row = row.to_dict()
    with np.errstate(divide='ignore'):
        row['score'] = int(min(1000, -100 * np.log10(row['wass_fdr'])))
    row['start'] = int(row['start'])
    row['end'] = int(row['end'])
    row['nreads_cntrl'] = int(row['nreads_cntrl'])
    row['nreads_treat'] = int(row['nreads_treat'])
    record = (
        '{chrom}\t{start:d}\t{end:d}\t'
        '{locus_id}\t{score:d}\t{strand}\t'
        '{wass_dist:.1f}\t{wass_dir:.1f}\t'
        '{wass_pval:.3g}\t{wass_fdr:.3g}\t'
        '{silhouette:.3f}\t{sil_pval:.3g}\t{sil_fdr:.3g}\t'
        '{nreads_cntrl:d}\t{nreads_treat:d}\n'
    )
    return record.format(**row)


def write_wass_test_output_bed(output_bed_fn, results):
    with open(output_bed_fn, 'w') as bed:
        for _, row in results.iterrows():
            record = format_bed_record(row)
            bed.write(record)