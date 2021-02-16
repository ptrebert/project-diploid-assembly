#!/usr/bin/env python3

import os
import io
import argparse

import pandas as pd


def parse_command_line():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--contig-table',
        '-ct',
        type=str,
        dest='contig_table'
    )
    parser.add_argument(
        '--fasta-folder',
        '-ff',
        type=str,
        dest='fasta_folder'
    )
    parser.add_argument(
        '--output-folder',
        '-of',
        type=str,
        dest='output_folder'
    )
    args = parser.parse_args()
    return args


def read_seqs_from_fasta(fasta_path, contigs):

    seq_buffer = io.StringIO()

    buffer = False
    with open(fasta_path, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                this_contig = line.strip().strip('>')
                if this_contig in contigs:
                    _ = seq_buffer.write('\n')
                    _ = seq_buffer.write(line)
                    buffer = True
                else:
                    buffer = False
                continue
            if buffer:
                _ = seq_buffer.write(line)
    return seq_buffer


def cache_fasta_paths(fasta_folder):

    cache = dict()
    for filename in os.listdir(fasta_folder):
        if not filename.endswith('.fasta'):
            continue
        sample, _, platform, _ = filename.split('_', 3)
        if platform == 'pbsq2-clr':
            tech = 'CLR'
        elif platform == 'pbsq2-ccs':
            tech = 'CCS'
        else:
            raise ValueError(filename)
        if 'h1-un' in filename:
            hap = 'H1'
        elif 'h2-un' in filename:
            hap = 'H2'
        else:
            raise ValueError(filename)
        cache[(sample, hap, tech)] = filename
    return cache


def main():

    args = parse_command_line()
    os.makedirs(args.output_folder, exist_ok=True)

    df = pd.read_csv(args.contig_table, sep='\t', header=0)

    fasta_cache = cache_fasta_paths(args.fasta_folder)

    for (sample, hap, tech), contigs in df.groupby(['sample', 'haplotype', 'platform']):
        try:
            fasta_file = fasta_cache[(sample, hap, tech)]
        except KeyError:
            print('skipping ', sample, hap, tech)
            continue
        contig_names = set(contigs['contig_id'])
        fasta_path = os.path.join(args.fasta_folder, fasta_file)
        contig_seqs = read_seqs_from_fasta(fasta_path, contig_names)

        outname = fasta_file.replace('.fasta', '.ctg3q29.fasta')
        outpath = os.path.join(args.output_folder, outname)
        if os.path.isfile(outpath):
            continue

        with open(outpath, 'w') as dump:
            _ = dump.write(contig_seqs.getvalue())
        print('done ', sample, hap, tech)
        
    return 0


if __name__ == '__main__':
    main()