#!/usr/bin/env python3

import os
import csv
import argparse as argp
import itertools as itt
import subprocess as sp
import functools
import multiprocessing as mp


JACCARD_CMD = 'bedtools jaccard -a {file_a} -b {file_b}'


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
    parser.add_argument(
        '--num-cpu',
        '-n',
        type=int,
        default=1,
        dest='num_cpu'
    )
    args = parser.parse_args()
    return args


@functools.lru_cache(128)
def parse_filename(file_path):

    fname = os.path.basename(file_path)
    sample, short_reads, long_reads, polisher, annotation, match_ratio, hap, ext = fname.split('.')
    platform = long_reads.split('_')[1].split('-')[1]
    assert platform in ['clr', 'ccs'], 'Unknown long read tech: {}'.format(long_reads)
    platform = 'HiFi' if platform == 'ccs' else 'CLR'
    assert hap in ['hap1-only', 'hap2-only'], 'Unknown haplotype: {}'.format(hap)
    assembly = 'HAP'
    haplotype = 1 if hap == 'hap1-only' else 2
    return sample, platform, assembly, haplotype


def collect_input_files(input_paths):

    input_files = []
    for ip in input_paths:
        if os.path.isdir(ip):
            files_under_path = [os.path.join(ip, f) for f in os.listdir(ip) if f.endswith('.bed')]
            input_files.extend(files_under_path)
        elif os.path.isfile(ip):
            input_files.append(ip)
        else:
            raise ValueError('Neither folder nor file: {}'.format(ip))
    assert input_files, 'No input files collected'
    return sorted(input_files)


def parse_jaccard_output(call_out):

    call_out = call_out.decode('utf-8')
    _, values = call_out.strip().split('\n')
    values = values.strip().split()
    header = [
        'intersect_bp',
        'union_bp',
        'jaccard',
        'intersect_count'
        ]
    converter = [int, int, float, int]
    parsed_output = dict((k, c(v)) for k, v, c in zip(header, values, converter))
    return parsed_output


def compute_jaccard(file_pair):

    file_a, file_b = file_pair

    sample_header = [
        'sample',
        'platform',
        'assembly',
        'haplotype'
    ]

    infos_a = parse_filename(file_a)
    infos_b = parse_filename(file_b)
    try:
        cmd_out = sp.check_output(
            JACCARD_CMD.format(**{
                'file_a': file_a,
                'file_b': file_b
            }),
            stderr=sp.STDOUT,
            timeout=60,
            shell=True,
            )
        
    except sp.CalledProcessError as spe:
        raise RuntimeError('Could not process {} / {}: {} / {}'.format(file_a, file_b, spe.returncode, str(spe)))

    infos_jaccard = parse_jaccard_output(cmd_out)
    row = dict((h + '_a', v) for h, v in zip(sample_header, infos_a))
    row.update(dict((h + '_b', v) for h, v in zip(sample_header, infos_b)))
    row.update(infos_jaccard)

    return row


def main():
    args = parse_command_line()

    input_files = collect_input_files(args.input)

    rows = []
    with mp.Pool(args.num_cpu) as pool:
        resit = pool.imap_unordered(compute_jaccard, itt.combinations(input_files, 2))
        [rows.append(res) for res in resit]

    header = sorted(rows[0].keys())

    with open(args.output, 'w') as dump:
        dictwriter = csv.DictWriter(dump, fieldnames=header, delimiter='\t')
        dictwriter.writeheader()
        dictwriter.writerows(rows)
    
    return 0


if __name__ == '__main__':
    main()