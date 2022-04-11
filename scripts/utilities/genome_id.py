#!/usr/bin/env python3

import sys
import argparse as argp
import pathlib as pl
import hashlib as hl
import csv
import gzip
import bz2
import lzma


MAGIC_NUMBERS = {
    bytes([0x42, 0x5A, 0x68]): ('bzip2', bz2.open, 'rt'),
    bytes([0x1F, 0x9D]): ('tar/zip (LZW)', None, None),
    bytes([0x1F, 0xA0]): ('tar/zip (LZH)', None, None),
    bytes([0x4C, 0x5A, 0x49, 0x50]): ('lzip', lzma.open, 'rt'),
    bytes([0x50, 0x4B, 0x03, 0x04]): ('zip', None, None),
    bytes([0x1F, 0x8B]): ('gzip (bgzip)', gzip.open, 'rt'),
    bytes([0x37, 0x7A, 0xBC, 0xAF, 0x27, 0x1C]): ('7-zip', None, None),
    bytes([0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00]): ('XZ/LZMA', lzma.open, 'rt'),
}

MAX_HEADER_LENGTH = max([len(k) for k, v in MAGIC_NUMBERS.items() if v[1] is not None])

NUCLEOTIDE_COMPLEMENT_MAP = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}

COMPLEMENT_TRANSLATION_TABLE = str.maketrans(NUCLEOTIDE_COMPLEMENT_MAP)


def parse_command_line():

    parser = argp.ArgumentParser(prog='genome_id.py', add_help=True)
    excl_group = parser.add_mutually_exclusive_group()
    excl_group.add_argument(
        '--describe',
        '--desc',
        '-d',
        action='store_true',
        default=False,
        dest='describe'
    )
    excl_group.add_argument(
        '--compare',
        '--comp',
        '-c',
        action='store_true',
        default=False,
        dest='compare'
    )
    parser.add_argument(
        '--genome',
        '-g',
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest='genome',
    )
    parser.add_argument(
        '--output',
        '-o',
        type=lambda x: pl.Path(x).resolve(),
        dest='output',
    )
    parser.add_argument(
        '--id-a',
        '-a',
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest='id_a',
    )
    parser.add_argument(
        '--id-b',
        '-b',
        type=lambda x: pl.Path(x).resolve(strict=True),
        dest='id_b',
    )
    parser.add_argument(
        '--lenient',
        '-l',
        action='store_true',
        default=False,
        dest='lenient'
    )
    parser.add_argument(
        '--report',
        '-r',
        action='store_true',
        default=False,
        dest='report'
    )
    args = parser.parse_args()
    if not args.compare:
        setattr(args, 'describe', True)

    return args



def determine_file_loader(genome_file):

    format_name = 'unknown'
    open_fun = None
    read_mode = None

    try:
        with open(genome_file, 'r') as fasta:
            for line in fasta:
                if not line.strip():
                    continue
                if line[0] == '>':
                    open_fun = open
                    read_mode = 'r'
                    format_name = 'plain/text'
                break
    except UnicodeDecodeError:
        pass

    if open_fun is None:
        with open(genome_file, 'rb') as compressed:
            header = compressed.read(MAX_HEADER_LENGTH)
            for magic_number, format_desc in MAGIC_NUMBERS.items():
                if header.startswith(magic_number):
                    format_name, open_fun, read_mode = format_desc
                    break

    if open_fun is None:
        raise RuntimeError(f'Genome file format is {format_name}, which is not supported')

    return open_fun, read_mode


def revcomp(sequence):
    return sequence.translate(COMPLEMENT_TRANSLATION_TABLE)[::-1]


def select_canonical_sequence(sequence):
    rc = revcomp(sequence)
    if rc < sequence:
        canonical = rc
    else:
        canonical = sequence
    return canonical


def describe_genome(genome_file):
    
    open_fun, read_mode = determine_file_loader(genome_file)

    genome_desc = []
    with open_fun(genome_file, read_mode) as fasta:
        seq_buffer = ''
        seq_name = None
        for line in fasta:
            if line.startswith('>'):
                if seq_name is not None:
                    seq_length = len(seq_buffer)
                    canonical_sequence = select_canonical_sequence(seq_buffer.upper())
                    seq_sha256 = hl.sha256(canonical_sequence.encode('ascii')).hexdigest()
                    genome_desc.append((seq_length, seq_name, seq_sha256))
                # in case comments are part of the FASTA header,
                # do not include them in the sequence name
                seq_name = line.strip()[1:].split()[0]
                seq_buffer = ''
                continue
            seq_buffer += line.strip()

    if seq_name is not None:
        seq_length = len(seq_buffer)
        canonical_sequence = select_canonical_sequence(seq_buffer.upper())
        seq_sha256 = hl.sha256(canonical_sequence.encode('ascii')).hexdigest()
        genome_desc.append((seq_length, seq_name, seq_sha256))

    if not genome_desc:
        raise ValueError(f'No sequence records read from genome file {genome_file} - is that a FASTA file?')

    return sorted(genome_desc, reverse=True)


def load_genome_id(genome_id_file):

    genome_id_strict = set()
    genome_id_lenient = set()
    with open(genome_id_file, 'r') as table:
        for line in table:
            seq_length, seq_name, seq_hash = line.strip().split('\t')
            genome_id_strict.add((seq_length, seq_name, seq_hash))
            genome_id_lenient.add((seq_length, seq_hash))
    return genome_id_strict, genome_id_lenient


def compare_genomes(genome_a, genome_b, lenient):

    a_strict, a_lenient = load_genome_id(genome_a)
    b_strict, b_lenient = load_genome_id(genome_b)

    a_only = a_strict - b_strict
    b_only = b_strict - a_strict
    if a_only or b_only:
        if lenient:
            a_only = a_lenient - b_lenient
            b_only = b_lenient - a_lenient
            if a_only or b_only:
                exit_code = 1
                a_only = '\n'.join(['\t'.join(t) for t in sorted(a_only, reverse=True)])
                b_only = '\n'.join(['\t'.join(t) for t in sorted(b_only, reverse=True)])
                report = '=== Genomes are not identical (lenient check) ===\n'
                report += f'Only in A ({genome_a.name}):\n'
                report += f'{a_only}\n'
                report += f'Only in B ({genome_b.name}):\n'
                report += f'{b_only}\n'
            else:
                report = 'Genomes have identical sequence content but different sequence names'
                exit_code = 0
        else:
            exit_code = 1
            a_only = '\n'.join(['\t'.join(t) for t in sorted(a_only, reverse=True)])
            b_only = '\n'.join(['\t'.join(t) for t in sorted(b_only, reverse=True)])
            report = '=== Genomes are not identical (strict check) ===\n'
            report += f'Only in A ({genome_a.name}):\n'
            report += f'{a_only}\n'
            report += f'Only in B ({genome_b.name}):\n'
            report += f'{b_only}\n'
    else:
        report = ''
        exit_code = 0
    return exit_code, report


def main():
    exit_code = 1
    args = parse_command_line()
    if args.describe:
        assert args.genome.is_file()
        genome_desc = describe_genome(args.genome)
        args.output.parent.mkdir(exist_ok=True, parents=True)
        with open(args.output, 'w') as table:
            writer = csv.writer(
                table,
                delimiter='\t',
                quoting=csv.QUOTE_NONE
            )
            writer.writerows(genome_desc)
        exit_code = 0
    elif args.compare:
        exit_code, report = compare_genomes(args.id_a, args.id_b, args.lenient)
        if args.report and report:
            sys.stderr.write(f'\n{report}\n')
    else:
        raise RuntimeError('Mode of operatio is not specified')

    return exit_code


if __name__ == '__main__':
    sys.exit(main())
