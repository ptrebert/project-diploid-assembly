#!/usr/bin/env python3

import os
import io
import collections as col
import argparse as argp
import multiprocessing as mp

import pandas as pd


DEFAULT_REFERENCE_SIZE = int(3.1e9)


def parse_command_line():
    parser = argp.ArgumentParser()

    parser.add_argument(
        '--assembly-fai',
        '-a',
        type=str,
        required=True,
        dest='assembly'
    )
    parser.add_argument(
        '--reference-fai',
        '-r',
        default=None,
        dest='reference'
    )
    parser.add_argument(
        '--output',
        '-o',
        type=str,
        required=True,
        dest='output'
    )
    args = parser.parse_args()
    return args


def load_reference_genome_size(ref_fai, default=DEFAULT_REFERENCE_SIZE):
    if isinstance(ref_fai, int):
        ref_size = ref_fai
        ref_name = 'custom'
    elif isinstance(ref_fai, float):
        ref_size = int(round(ref_fai, 0))
        ref_name = 'custom'
    elif os.path.isfile(ref_fai):
        df = pd.read_csv(ref_fai, sep='\t', usecols=[0,1], names=['chrom', 'size'], header=None)
        ref_name = os.path.basename(ref_fai).rsplit('.', 2)[0]
        ref_size = int(df['size'].sum())
    else:
        ref_size = default
        ref_name = 'custom'
    assert ref_size > 0, 'Reference genome size 0: {} / default: {}'.format(ref_fai, default)
    ref_size_gbp = round(ref_size / 1e9, 2)
    return ref_size, ref_size_gbp, ref_name


def compute_assembly_stats(contigs, stats, ref_size_bp):

    stats.append(('assembly_contigs', contigs.shape[0]))
    assm_size_bp = int(contigs['size'].sum())
    assm_size_gbp = round(assm_size_bp / 1e9, 2)
    stats.append(('assembly_size_bp', assm_size_bp))
    stats.append(('assembly_size_Gbp', assm_size_gbp))

    mean_ctg_length = int(contigs['size'].mean())
    median_ctg_length = int(contigs['size'].median())
    max_ctg_length = int(contigs['size'].max())
    min_ctg_length = int(contigs['size'].min())

    stats.extend([
        ('assembly_contig_length_min', min_ctg_length),
        ('assembly_contig_length_median', median_ctg_length),
        ('assembly_contig_length_mean', mean_ctg_length),
        ('assembly_contig_length_max', max_ctg_length)
    ])

    n50_threshold = assm_size_bp // 2
    n90_threshold = int(round(assm_size_bp * 0.9, 0))

    contigs.sort_values('size', inplace=True, axis=0, ascending=False)
    contigs['NX'] = contigs['size'].cumsum()

    n50_value = contigs.loc[contigs['NX'] >= n50_threshold, 'size'].max()
    n90_value = contigs.loc[contigs['NX'] >= n90_threshold, 'size'].max()

    stats.append(('assembly_contig_length_N50', n50_value))
    stats.append(('assembly_contig_length_N90', n90_value))

    reference_cov = round(assm_size_bp / ref_size_bp, 2)
    stats.append(('assembly_reference_coverage', reference_cov))

    return stats


def main():
    args = parse_command_line()

    ref_size_bp, ref_size_gbp, ref_name = load_reference_genome_size(args.reference)

    stats = [
        ('reference_genome', ref_name),
        ('reference_size_bp', ref_size_bp),
        ('reference_size_Gbp', ref_size_gbp),
    ]

    contigs = pd.read_csv(
        args.assembly,
        sep='\t',
        usecols=[0,1],
        names=['contig', 'size'],
        header=None
    )
    stats.append(('assembly_name', os.path.basename(args.assembly).rsplit('.', 2)[0]))

    stats = compute_assembly_stats(contigs, stats, ref_size_bp)

    stats = pd.DataFrame.from_records(stats, columns=['statistic', 'value'])

    stats.to_csv(args.output, sep='\t', index=False, header=False)

    return 0


if __name__ == '__main__':
    main()
