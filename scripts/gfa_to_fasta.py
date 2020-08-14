#!/usr/bin/env python3

import os
import io
import argparse as argp
import multiprocessing as mp

import pandas as pd


def parse_command_line():
    parser = argp.ArgumentParser()

    parser.add_argument(
        '--gfa',
        '-g',
        type=str,
        required=True,
        dest='gfa'
    )
    parser.add_argument(
        '--n-cpus',
        '-n',
        type=int,
        default=1,
        dest='n_cpus'
    )
    parser.add_argument(
        '--out-fasta',
        '-of',
        type=str,
        required=True,
        dest='out_fasta'
    )
    parser.add_argument(
        '--out-map',
        '-om',
        type=str,
        required=True,
        dest='out_map'
    )
    parser.add_argument(
        '--out-stats',
        '-os',
        type=str,
        required=True,
        dest='out_stats'
    )
    args = parser.parse_args()
    return args


def parse_gfa_segment(line):
    """
    Primarily intended for hifiasm segment lines
    """
    columns = line.strip().split()
    contig = columns[1]
    seq = columns[2]
    seqlen = int(columns[-2].split(':')[-1])
    assert seqlen == len(seq), 'length mismatch: {} / {} / {}'.format(contig, seqlen, len(seq))
    stats_counter = col.Counter(seq)
    stats_counter['LEN'] = seqlen
    stats = dict(stats_counter)
    stats['ID'] = contig
    
    buffer = io.StringIO()
    _ = buffer.write('>{}\n'.format(contig))
    bases_buffered = 0
    for i in range(0, seqlen // 120 + 1):
        buffered = buffer.write(seq[i*120:i*120+120] + '\n')
        bases_buffered += (buffered - 1)
    assert bases_buffered == seqlen, 'Dropped bases for {}: {} / {}'.format(contig, seqlen, bases_buffered)
    return buffer, stats_counter


def parse_gfa_assembly(line):
    columns = line.strip().split()
    contig = columns[1]
    start = int(columns[2])
    orient = columns[3]
    end = start + int(columns[6])
    name = columns[4]
    return contig, start, end, name, orient


def parse_gfa_line(line):
    if line.startswith('S'):
        return parse_gfa_segment(line)
    elif line.startswith('A'):
        return parse_gfa_assembly(line)
    else:
        raise ValueError('Unexpected GFA line: {} .....'.format(line.strip()[:50]))


def convert_gfa_to_fasta(gfa_file, threads):
    """
    Primarily intended to convert hifiasm gfa to FASTA
    """

    fasta_buffer = io.StringIO()
    read_to_contigs = []
    contig_stats = []

    with open(gfa_file, 'r') as gfa:

        with mp.Pool(threads) as pool:

            res_iter = pool.imap_unordered(parse_gfa_line, gfa)
            for res in res_iter:
                if isinstance(res, tuple) and len(res) == 5:
                    read_to_contigs.append(res)
                elif isinstance(res, tuple) and len(res) == 2:
                    fasta_buffer.write(res[0].getvalue())
                    contig_stats.append(res[1])
                else:
                    raise ValueError('Unexpected result gfa-to-fasta: {}'.format(res))

    return fasta_buffer, read_to_contigs, contig_stats


def main():
    args = parse_command_line()

    fasta_buffer, rc_map, stats = convert_gfa_to_fasta(args.gfa, args.ncpus)

    with open(args.out_fasta, 'w') as dump:
        _ = dump.write(fasta_buffer.getvalue())
        
    df_map = pd.DataFrame.from_records(
        rc_map,
        columns=['contig', 'start', 'end', 'read', 'orientation'])
    df_map.sort_values(['contig', 'start'], inplace=True)
    df_map.to_csv(args.out_map, sep='\t', header=True, index=False)

    df_stats = pd.DataFrame.from_records(
        stats,
        columns=['ID', 'LEN', 'A', 'C', 'G', 'T']
    )
    df_stats.sort_values(['LEN', 'ID'], inplace=True)
    df_stats.to_csv(args.out_stats, sep='\t', header=True, index=False)

    return 0


if __name__ == '__main__':
    main()
