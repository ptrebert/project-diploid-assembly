#!/usr/bin/env python3

import os
import argparse

import pandas as pd


def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--bed-folder',
        '-b',
        dest='bed_folder',
        type=str,
    )
    parser.add_argument(
        '--chrom-sizes',
        '-c',
        dest='chrom_sizes',
        type=str
    )
    parser.add_argument(
        '--sample-table',
        '-s',
        type=str,
        dest='sample_table'
    )
    parser.add_argument(
        '--output',
        '-o',
        dest='output',
        type=str
    )
    args = parser.parse_args()
    return args


def load_tabular_data(file_path, header='infer', names=None):

    if names is not None:
        df = pd.read_csv(file_path, sep='\t', names=names, header=0)
    else:
        df = pd.read_csv(file_path, sep='\t', header=header)
    return df


def load_chromosome_table(file_path):

    chrom_table = load_tabular_data(file_path, header=None)
    chrom_table = chrom_table[chrom_table.columns[:2]]
    chrom_table.columns = ['chrom', 'size']
    chrom_table['size'] = chrom_table['size'].astype('int32')

    return chrom_table


def load_sample_table(file_path):

    sample_table = load_tabular_data(file_path)

    # ? TODO ?
    # drop "skip" samples...

    return sample_table


def extract_sample_info(file_name):

    sample = file_name.split('_')[0]
    if 'pbsq2-clr' in file_name:
        tech = 'CLR'
    elif 'pbsq2-ccs' in file_name:
        tech = 'HiFi'
    else:
        raise ValueError('Cannot determine read type: {}'.format(file_name))

    if 'h1-un' in file_name:
        hap = 'H10'
    elif 'h2-un' in file_name:
        hap = 'H20'
    else:
        raise ValueError('Cannot determine haplotype: {}'.format(file_name))
    
    return sample, tech, hap


def load_contig_alignments(bed_folder, sample_table):

    if 'sample' in sample_table:
        sample_column = 'sample'
    elif 'individual' in sample_table:
        sample_column = 'individual'
    else:
        raise ValueError('Cannot identify sample column in table: {}'.format(sample_table.head()))

    orientation_map = {
        '+': 1,
        '-': -1,
        '.': 0
    }

    store_values = []

    for root, dirs, files in os.walk(bed_folder, followlinks=False):
        files = [f for f in files if f.endswith('.bed')]
        for f in files:
            sample, technology, haplotype = extract_sample_info(f)
            super_pop = sample_table.loc[sample_table[sample_column] == sample, 'super_population']
            assert len(super_pop) == 1, 'Multi super pop: {}'.format(super_pop)
            super_pop = super_pop.values[0]
            store_key = os.path.join(super_pop, sample, technology, haplotype)

            aln = load_tabular_data(
                os.path.join(root, f),
                0,
                ['chrom', 'start', 'end', 'contig_name', 'mapq', 'orientation']
            )
            aln['orientation'] = aln['orientation'].apply(lambda x: orientation_map[x])
            aln['orientation'] = aln['orientation'].astype('int8')
            aln['start'] = aln['start'].astype('int32')
            aln['end'] = aln['end'].astype('int32')
            aln['mapq'] = aln['mapq'].astype('int8')
            store_values.append((store_key, aln))

    return store_values
            

def main():
    args = parse_command_line()

    chrom_table = load_chromosome_table(args.chrom_sizes)
    sample_table = load_sample_table(args.sample_table)

    contig_alignments = load_contig_alignments(
        args.bed_folder,
        sample_table
    )

    with pd.HDFStore(args.output, mode='w', complevel=9, complib='blosc') as hdf:
        hdf.put(os.path.join('metadata', 'chrom'), chrom_table, format='fixed')
        hdf.put(os.path.join('metadata', 'sample'), sample_table, format='fixed')

        for k, v in contig_alignments:
            hdf.put(k, v, format='fixed')

    return


if __name__ == '__main__':
    main()
