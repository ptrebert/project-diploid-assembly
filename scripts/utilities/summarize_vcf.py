#!/usr/bin/env python3

import os
import argparse as argp

import pandas as pd
import pysam


GT_MAP = {
    (0, 0): 'HOM',
    (1, 1): 'HOM',
    (0, 1): 'HET',
    (1, 0): 'HET'
}


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        '--input',
        '-i',
        type=str,
        dest='input',
        nargs='+'
        )
    parser.add_argument(
        '--output',
        '-o',
        type=str,
        dest='output',
        required=True
    )
    args = parser.parse_args()
    return args


def parse_filename(file_path):

    fname = os.path.basename(file_path)
    short_reads, long_read_assm = fname.split('_map-to_')
    sample = short_reads.split('_')[0]

    hap = long_read_assm.split('.')[1]
    if hap == 'h1-un':
        hap = 10
    elif hap == 'h2-un':
        hap = 20
    else:
        raise ValueError('Unexpected haplotype: {}'.format(fname))

    platform = long_read_assm.split('_')[1].split('-')[1]
    assert platform in ['clr', 'ccs'], 'Unknown long read tech: {}'.format(long_reads)
    platform = 'HiFi' if platform == 'ccs' else 'CLR'

    return sample, platform, hap


def main():
    args = parse_command_line()

    out_mode = 'w'

    for vcf_file in args.input:
        rows = []
        index = []
        sample, platform, hap = parse_filename(vcf_file)

        with pysam.VariantFile(vcf_file, 'r') as vcf:
            for record in vcf:
                assert record.chrom == record.contig, 'Sequence mismatch: {}'.format(record)
                v = {
                    'sequence': record.chrom,
                    'pos': record.pos,
                    'start': record.start,
                    'stop': record.stop,
                    'qual': int(round(record.qual, 0)),
                    'ref_allele': record.ref,
                    'alt_allele': record.alts[0],
                    'depth': record.info['DP'],
                    'region_length': record.rlen,
                    'variant_length': record.info['LEN'][0]
                }
                var_type = record.info['TYPE'][0].upper()
                if var_type == 'SNP':
                    var_type = 'SNV'
                genotype = record.samples[sample]['GT']
                gt = GT_MAP[genotype]
                index.append((sample, platform, hap, var_type, gt))
                rows.append(v)
    
        df = pd.DataFrame.from_records(
            rows,
            index=pd.MultiIndex.from_tuples(
                index,
                names=['sample', 'platform', 'hap', 'var_type', 'genotype']
            )
        )
        store_key = os.path.join(sample, platform, 'HAP' + str(hap))
        df.to_hdf(args.output, store_key, mode=out_mode, format='fixed', complevel=9)
        out_mode = 'a'

    return 0


if __name__ == '__main__':
    main()