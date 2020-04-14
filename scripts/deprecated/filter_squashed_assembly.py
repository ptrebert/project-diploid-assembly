#!/usr/bin/env python3

import os as os
import sys as sys
import logging as log
import io as io
import gzip as gz
import traceback as trb
import argparse as argp
import pickle as pck
import time as ti
import re as re
import functools as fnt
import collections as col
import multiprocessing as mp


import numpy as np
import pysam as pysam


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog="filter_squashed_assembly.py", description=__doc__)
    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        default=False,
        dest="debug",
        help="Print status and error messages to STDOUT. Otherwise, only"
        " errors/warnings will be reported to STDERR.",
    )
    parser.add_argument(
        "--input-fasta",
        "-if",
        type=str,
        required=True,
        dest="input",
        help="Full path FASTA file (uncompressed).",
    )
    parser.add_argument(
        "--output-fasta",
        "-of",
        type=str,
        required=True,
        dest="output",
        help="Full path to output FASTA file.",
    )
    parser.add_argument(
        "--output-metrics",
        "-om",
        type=str,
        required=True,
        dest="metrics",
        help="Full path to metrics (statistics) file."
    )
    parser.add_argument(
        "--output-regions",
        "-or",
        type=str,
        required=False,
        default="",
        dest="regions",
        help="Full path to regions file. Contigs are split into"
             " --min-size'd regions (0-based, half-open)."
    )
    parser.add_argument(
        "--min-size",
        "-ms",
        type=int,
        default=100000,
        dest='minsize',
        help="Minimum size of contigs to be kept in output."
             " Default: 100000"
    )
    return parser.parse_args()


def read_fasta_input(fpath, logger):
    """
    :param fpath:
    :param logger:
    :return:
    """
    logger.debug('Reading input FASTA from path {}'.format(fpath))

    contig_seqs = dict()
    cluster_to_contig = dict()
    contig_sizes = []

    contig_name = None
    contig_size = None
    contig_splits = None
    seq_buffer = None
    nucleotides = col.Counter()

    with open(fpath, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                if seq_buffer is not None:
                    if contig_splits is not None:
                        assert contig_size > 0, 'Did not record contig size'
                        contig_sizes.append((contig_size, contig_name))

                    contig_seqs[contig_name] = seq_buffer.getvalue()
                    nucleotides.update(col.Counter(seq_buffer.getvalue().upper()))
                contig_size = 0
                contig_splits = None
                seq_buffer = io.StringIO()
                try:
                    contig_name, contig_size = line.strip().strip('>').split()
                    # format: wtdbg2 output
                    contig_size = int(contig_size.strip('len='))
                    contig_sizes.append((contig_size, contig_name))
                except ValueError:
                    logger.debug('Detected non-wtdbg2 format, assuming SaarClust-v0 assembly')
                    try:
                        contig_name, contig_splits = line.strip().strip('>').split('_')
                    except ValueError:
                        logger.debug('Nope... assuming SaarClust-v1 assembly')
                        contig_name = line.strip().strip('>')
                        if '.' in contig_name:
                            contig_name, contig_splits = contig_name.split('.', 1)
                            contig_splits = contig_splits.replace('.', '-')
                        else:
                            contig_splits = 'contig-to-cluster-unknown'
                        cluster_to_contig[contig_name] = contig_splits
                        contig_size = 0
                    else:
                        contig_splits = contig_splits.split('.')
                        contig_splits = sorted(contig_splits, key=lambda x: int(x.strip('ctg')))
                        cluster_to_contig[contig_name] = contig_splits
                        contig_size = 0
            else:
                contig_size += len(line.strip())
                seq_buffer.write(line)

    contig_seqs[contig_name] = seq_buffer.getvalue()
    nucleotides.update(col.Counter(seq_buffer.getvalue().upper()))
    if contig_splits is not None:
        contig_sizes.append((contig_size, contig_name))

    logger.debug('Read {} contigs from FASTA input'.format(len(contig_sizes)))

    contig_sizes = sorted(contig_sizes, reverse=True)

    return contig_sizes, contig_seqs, nucleotides, cluster_to_contig


def write_fasta_output(fpath, minsize, contig_sizes, contig_seqs, logger):
    """
    :param fpath:
    :param minsize:
    :param contig_sizes:
    :param contig_seqs:
    :param logger:
    :return:
    """
    logger.debug('Creating output path...')
    os.makedirs(os.path.abspath(os.path.dirname(fpath)), exist_ok=True)

    with open(fpath, 'w') as fasta:
        for length, name in contig_sizes:
            if length < minsize:
                logger.debug('Hit size limit: {} / {}'.format(name, length))
                break
            _ = fasta.write('>{}\n'.format(name))
            _ = fasta.write(contig_seqs[name])
    logger.debug('Output FASTA written to {}'.format(fpath))
    return


def main(logger, cargs):
    """
    :param logger:
    :param cargs:
    :return:
    """
    logger.debug("Starting computations")

    contig_sizes, contig_seqs, nucleotides, cluster_to_contig = read_fasta_input(cargs.input, logger)

    metrics_out = io.StringIO()
    metrics_out.write('input_fasta\t{}\n'.format(cargs.input))

    logger.debug('Writing nucleotide counts')

    total_nucleotides = nucleotides['A'] + nucleotides['C'] + nucleotides['G'] + nucleotides['T'] + nucleotides['N']
    for nuc in ['A', 'C', 'G', 'T']:
        metrics_out.write('nucleotide_fraction_{}\t{}\n'.format(nuc, round(nucleotides[nuc] / total_nucleotides, 3)))

    if len(cluster_to_contig) > 0:
        try:
            sorted_clusters = sorted(cluster_to_contig.keys(), key=lambda x: int(x.strip('cluster')))
        except ValueError:
            sorted_clusters = sorted(cluster_to_contig.keys(), key=lambda x: 23 if x == 'chrX' else int(x.strip('chr')))
        for cluster in sorted_clusters:
            cluster_contigs = cluster_to_contig[cluster]
            if 'unknown' in cluster_contigs or '-' in cluster_contigs:
                pass
            else:
                cluster_contigs = ','.join(cluster_to_contig[cluster])
            metrics_out.write('{}\t{}\n'.format(cluster, cluster_contigs))

    logger.debug('Computing FASTA input statistics')

    total_seq_len = sum([l for l, n in contig_sizes])
    assert total_seq_len == total_nucleotides, \
        'Seq. len. / nucleotide mismatch: {} v {}'.format(total_seq_len, total_nucleotides)
    metrics_out.write('sequence_bp_total\t{}\n'.format(total_seq_len))

    seq_len_at_minsize = sum([l for l, n in contig_sizes if l >= cargs.minsize])
    metrics_out.write('sequence_bp_geq_{}\t{}\n'.format(cargs.minsize, seq_len_at_minsize))
    metrics_out.write('sequence_fraction_geq_{}\t{}\n'.format(cargs.minsize, round(seq_len_at_minsize / total_seq_len, 3)))

    half_minsize = int(0.5 * cargs.minsize)
    seq_len_at_half_minsize = sum([l for l, n in contig_sizes if l >= half_minsize])
    metrics_out.write('sequence_bp_geq_{}\t{}\n'.format(half_minsize, seq_len_at_half_minsize))
    metrics_out.write('sequence_fraction_geq_{}\t{}\n'.format(half_minsize, round(seq_len_at_half_minsize / total_seq_len, 3)))

    double_minsize = 2 * cargs.minsize
    seq_len_at_double_minsize = sum([l for l, n in contig_sizes if l >= double_minsize])
    metrics_out.write('sequence_bp_geq_{}\t{}\n'.format(double_minsize, seq_len_at_double_minsize))
    metrics_out.write('sequence_fraction_geq_{}\t{}\n'.format(double_minsize, round(seq_len_at_double_minsize / total_seq_len, 3)))

    for length, name in contig_sizes:
        metrics_out.write('{}\t{}\n'.format(name, length))

    write_fasta_output(cargs.output, cargs.minsize, contig_sizes, contig_seqs, logger)

    logger.debug('Writing metrics file...')
    os.makedirs(os.path.abspath(os.path.dirname(cargs.metrics)), exist_ok=True)

    with open(cargs.metrics, 'w') as metrics:
        _ = metrics.write(metrics_out.getvalue())

    if cargs.regions:
        logger.debug('Preparing region splits')
        out_regions = []
        for length, name in contig_sizes:
            blank_region = name + ':{}-{}'
            splits = length // cargs.minsize
            remainder = length % cargs.minsize
            current_start = 0
            for split_num in range(splits):
                if split_num == (splits - 1):
                    if remainder < (cargs.minsize / 2):
                        out_regions.append(blank_region.format(current_start, length))
                    else:
                        out_regions.append(blank_region.format(current_start, current_start + cargs.minsize))
                        current_start += cargs.minsize
                        out_regions.append(blank_region.format(current_start, length))
                else:
                    out_regions.append(blank_region.format(current_start, current_start + cargs.minsize))
                    current_start += cargs.minsize
        logger.debug('Writing regions file...')
        os.makedirs(os.path.abspath(os.path.dirname(cargs.regions)), exist_ok=True)
        with open(cargs.regions, 'w') as regions:
            _ = regions.write('\n'.join(out_regions))

    return


if __name__ == "__main__":
    logger = None
    rc = 0
    try:
        log_msg_format = "%(asctime)s | %(levelname)s | %(message)s"
        cargs = parse_command_line()
        if cargs.debug:
            log.basicConfig(stream=sys.stdout, level=log.DEBUG, format=log_msg_format)
        else:
            log.basicConfig(stream=sys.stderr, level=log.WARNING, format=log_msg_format)
        logger = log.getLogger()
        logger.debug("Logger initiated")
        main(logger, cargs)
        logger.debug("Run completed - exit")
        log.shutdown()
    except Exception as exc:
        rc = 1
        if logger is not None:
            logger.error("Unrecoverable error: {}".format(str(exc)))
            logger.debug("=== TRACEBACK ===\n\n")
            buf = io.StringIO()
            trb.print_exc(file=buf)
            logger.error(buf.getvalue())
            logger.debug("Exit\n")
            log.shutdown()
        else:
            trb.print_exc()
    finally:
        sys.exit(rc)
