#!/usr/bin/env python3

import os as os
import sys as sys
import shutil as sh
import logging as log
import io as io
import csv as csv
import traceback as trb
import argparse as argp
import pickle as pck
import time as ti
import collections as col
import multiprocessing as mp


import numpy as np
import pysam as pysam
import dnaio as dnaio
import xopen as xopen


Read = col.namedtuple('Read', 'record_num name sequence')


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog="collect_read_stats.py", description=__doc__)
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
        "--input-files",
        "-if",
        type=str,
        required=True,
        dest="input",
        nargs='+',
        help="Full path FASTQ/FASTA read files, BAM alignments or 'seqtk comp' output "
            "(text-based files can be compressed).",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        dest="output",
        help="Full path to output file (Python pickle format). Directories will "
             "be created if they do not exist.",
    )
    parser.add_argument(
        "--summary-output",
        "-so",
        type=str,
        dest='summary_output',
        default=None,
        help="Full path to short summary output file listing the most important statistics. "
             "Requires specifying the genome size for coverage computation. Default: none"
    )
    parser.add_argument(
        "--buffer-size",
        "-bs",
        type=int,
        default=1000000,
        dest="buffer_size",
        help="Store this many G+C values before updating the binned statistics. Default: 1 000 000"
    )
    mutgrp = parser.add_mutually_exclusive_group()
    mutgrp.add_argument(
        "--genome-size-file",
        "-gsf",
        type=str,
        dest="size_file",
        default=None,
        help="Read genome size from this file (assumed: FASTA index), sum up all values in column 2."
    )
    mutgrp.add_argument(
        "--genome-size",
        "-gs",
        type=int,
        dest="genome_size",
        default=0,
        help="specify genome size (in bp). Default: 0"
    )
    parser.add_argument(
        "--num-cpu",
        "-n",
        type=int,
        default=1,
        dest="numcpu",
        help="Specify number of CPU cores to use for parallel computation. "
             "Note that there is no sanity check for actually free/available cores. "
             "Default: 1",
    )
    parser.add_argument(
        "--copy-stats-dump",
        "-cpd",
        type=str,
        nargs='*',
        default=None,
        dest="copy_dump",
        help="Path to file(s) containing the same data set statistics derived from a different "
             "input data format; simply copy instead of reprocessing."
    )
    parser.add_argument(
        "--copy-summary",
        "-cps",
        type=str,
        nargs='*',
        default=None,
        dest="copy_summary",
        help="Path to file(s) containing the same data set summary derived from a different "
             "input data format; simply copy instead of reprocessing."
    )
    return parser.parse_args()


def compute_read_statistics(read):
    """
    :param read:
    :return:
    """
    length = len(read.sequence)
    bases = col.Counter(read.sequence)
    pct_gc = round((bases['G'] + bases['C']) / length, 3)
    return length, bases, pct_gc


def extract_read_statistics(read):
    """
    """
    length = int(read.sequence['sequence_length'])
    bases = {
        'A': int(read.sequence['num_A']),
        'C': int(read.sequence['num_C']),
        'G': int(read.sequence['num_G']),
        'T': int(read.sequence['num_T']),
        'N': int(read.sequence['num_N'])
    }
    pct_gc = round((bases['G'] + bases['C']) / length, 3)
    return length, bases, pct_gc


def read_sequence_length_file(fpath):
    """
    :param fpath: FASTA index (.fai file)
    :return:
    """
    total_length = 0
    with open(fpath, 'r') as faidx:
        for line in faidx:
            if line.startswith('#'):
                continue
            seq_len = line.split()[1]
            try:
                seq_len = int(seq_len)
                total_length += seq_len
            except (ValueError, TypeError):
                raise ValueError('Malformed (non-integer) sequence length '
                                 'entry in file {}: {}'.format(fpath, line.strip()))
    return total_length


def read_text_sequence_records(fpath):
    """
    :param fpath:
    :return:
    """
    with dnaio.open(fpath) as fastx:
        for idx, record in enumerate(fastx, start=1):
            rd = Read(record_num=idx,
                      name=record.name,
                      sequence=record.sequence
                      )
            yield rd
    return


def read_seqtk_records(fpath):
    """
    :param fpath:
    :return:
    """
    header = [
        'sequence_name',
        'sequence_length',
        'num_A',
        'num_C',
        'num_G',
        'num_T',
        'num_RYKMSW',
        'num_BDHV',
        'num_N',
        'num_CpG',
        'num_transversions',
        'num_transitions',
        'num_CpG_transitions'
    ]
    with xopen.xopen(fpath, 'rt', threads=2) as table:
        reader = csv.DictReader(
            table,
            fieldnames=header,
            delimiter='\t'
        )
        for idx, row in enumerate(reader, start=1):
            rd = Read(
                record_num=idx,
                name=row['sequence_name'],
                sequence=row
            )
            yield rd
    return


def read_bam_sequence_records(fpath):
    """
    :param fpath:
    :return:
    """
    with pysam.AlignmentFile(fpath, 'rb', check_sq=False) as bam:
        for idx, record in enumerate(bam, start=1):
            rd = Read(record_num=idx,
                      name=record.query_name,
                      sequence=record.query_sequence
                      )
            yield rd
    return


def assemble_file_processors(input_files, logger):
    """
    """
    file_paths, chunk_readers, read_processors = [], [], []
    for input_file in input_files:
        file_paths.append(input_file)
        if any([x in input_file.lower() for x in ['.fastq', '.fq', '.fasta', '.fa']]):
            logger.debug('Adding FASTQ/FASTA processor for file: {}'.format(os.path.basename(input_file)))
            chunk_readers.append(read_text_sequence_records)
            read_processors.append(compute_read_statistics)
        elif '.bam' in input_file.lower():
            logger.debug('Adding BAM processors for file {}'.format(os.path.basename(input_file)))
            chunk_readers.append(read_bam_sequence_records)
            read_processors.append(compute_read_statistics)
        else:
            # assume seqtk comp output
            logger.debug('Adding seqtk processor for file {}'.format(os.path.basename(input_file)))
            chunk_readers.append(read_seqtk_records)
            read_processors.append(extract_read_statistics)
    return file_paths, chunk_readers, read_processors


def get_total_genome_size(fpath, gsize):
    """
    :param fpath:
    :param gsize:
    :return:
    """
    if (fpath is None and gsize == 0) or (fpath is not None and gsize > 0):
        raise ValueError('Summary output requires specifying genome size (either via file or via numeric value)')
    total_length = gsize
    if fpath is not None:
        total_length = read_sequence_length_file(fpath)
    return total_length


def determine_cov_steps(average_read_length, max_read_length):
    """
    Depending on the read length statistic, determine which coverage
    steps should be included in the textual summary output
    """
    if average_read_length < 1000:
        # ~ Illumina short
        cov_steps = list(range(0, max_read_length // 25 * 25 + 25, 25))
    elif average_read_length < 20000:
        # ~ HiFi reads
        cov_steps = [0, 1000]
        if max_read_length < 100000:
            cov_steps += list(range(5000, 55000, 5000))
        elif 100000 < max_read_length < 500000:
            cov_steps += [5000, 10000]
            cov_steps += list(range(20000, 220000, 20000))
        else:
            cov_steps += [5000, 10000, 25000]
            cov_steps += list(range(50000, 550000, 25000))
    elif average_read_length < 50000 and max_read_length < 200000:
        cov_steps = [0, 1000, 5000, 10000]
        cov_steps += list(range(20000, 110000, 10000))
    elif average_read_length < 50000 and max_read_length > 200000:
        cov_steps = [0, 1000, 5000, 10000]
        cov_steps += list(range(25000, 525000, 25000))
    else:
        cov_steps = [0, 1000, 5000, 20000]
        cov_steps += list(range(50000, int(1e6) + 50000, 50000))
    return cov_steps


def prepare_summary_statistics(length_stat, genome_size, num_reads):
    """
    :param length_stat:
    :param genome_size:
    :param num_reads:
    :return:
    """
    unpacked = np.array(
        [(l, length_stat[l]) for l in sorted(length_stat.keys(), reverse=True)],
        dtype=np.int32
    ).transpose()
    # 0: read length
    # 1: count/abundance

    unpacked = np.vstack((unpacked, unpacked[0, ] * unpacked[1, ]))
    # 2: coverage per length

    # cumulative read abundance
    cum_temp = unpacked[1, :].cumsum()
    temp_selector = np.asarray((cum_temp / cum_temp.max()) <= 0.5)
    median_read_length = unpacked[0, temp_selector].min()

    # cumulative coverage - used again below
    cum_temp = unpacked[2, :].cumsum()
    total_seq_bp = cum_temp.max()
    total_seq_giga_bp = round(total_seq_bp / 1e9, 2)
    temp_selector = np.asarray((cum_temp / cum_temp.max()) <= 0.5)
    n50_read_length = unpacked[0, temp_selector].min()
    max_read_length = unpacked[0, ].max()  # needed for coverage steps

    summary_stats = [
        ('num_reads', num_reads),
        ('total_seq_bp', total_seq_bp),
        ('total_seq_Gbp', total_seq_giga_bp),
        ('genome_size_bp', genome_size),
        ('genome_size_Gbp', round(genome_size / 1e9, 2)),
        ('read_length_min', unpacked[0, :].min()),
        ('read_length_median', median_read_length)
    ]

    # compute average read length via (abundance) weighted average
    avg_rlen = np.average(
        unpacked[0, ],
        weights=unpacked[1, ]
    )

    summary_stats.append(('read_length_mean', int(avg_rlen)))
    summary_stats.append(('read_length_max', max_read_length))
    summary_stats.append(('read_length_N50', n50_read_length))

    total_cov = round(total_seq_bp / genome_size, 2)
    summary_stats.append(('cov_total', total_cov))

    cov_steps = determine_cov_steps(avg_rlen, max_read_length)
    for s in cov_steps:
        temp_selector = np.asarray(unpacked[0, :] >= s)
        try:
            temp_cov = cum_temp[temp_selector].max()
            temp_cov = round(temp_cov / genome_size, 2)
        except ValueError:
            temp_cov = 0
        summary_stats.append(('cov_geq_{}'.format(s), temp_cov))

    return summary_stats


def compute_dataset_statistics(cargs, logger):
    """
    :param cargs:
    :param logger:
    :return:
    """
    file_paths, chunk_readers, read_processors = assemble_file_processors(cargs.input, logger)

    len_stats = col.Counter()
    nuc_stats = col.Counter()
    gc_stats = col.Counter()
    read_count = 0

    gc_bins = np.arange(0, 1.01, 0.01)
    # since numpy.digitize is not called as right-inclusive
    # adjust boundary to include 1
    gc_bins[-1] += 0.001

    gc_buffer = np.zeros(cargs.buffer_size, dtype=np.float16)
    gc_idx = 0

    with mp.Pool(cargs.numcpu) as pool:
        logger.debug('Initialized worker pool')
        for input_file, input_reader, read_processor in zip(file_paths, chunk_readers, read_processors):
            logger.debug('Processing file {}'.format(input_file))
            resit = pool.imap(read_processor, input_reader(input_file))
            for l, n, c in resit:
                read_count += 1
                len_stats[l] += 1
                nuc_stats.update(n)
                gc_buffer[gc_idx] = c
                gc_idx += 1
                if gc_idx == cargs.buffer_size:
                    logger.debug('Buffer limit reached - binning GC values')
                    gc_stats.update(col.Counter(np.digitize(gc_buffer, gc_bins, right=False)))
                    gc_buffer = np.zeros(cargs.buffer_size, dtype=np.float16)
                    gc_idx = 0
                if read_count % int(1e6) == 0:
                    logger.debug('Processed {} reads...'.format(read_count))
            logger.debug('All reads processed (total: {})'.format(read_count))

    if gc_idx != 0:
        gc_stats.update(col.Counter(np.digitize(gc_buffer[:gc_idx], gc_bins, right=False)))

    summary_stats = None
    total_genome_size = -1

    # the summary output could have been copied from a different source
    if cargs.summary_output is not None:
        if os.path.isfile(cargs.summary_output):
            logger.debug('Loading summary statistics from file: {}'.format(cargs.summary_output))
            summary_stats = {}
            with open(cargs.summary_output, 'r') as summary:
                for line in summary:
                    key, value = line.strip().split('\t')
                    if any([x in key for x in ['cov', 'Gbp']]):
                        value = float(value)
                    else:
                        value = int(value)
                    summary_stats[key] = value
        else:
            logger.debug('User requested summary output...')
            total_genome_size = get_total_genome_size(cargs.size_file, cargs.genome_size)
            gs_gbp = round(total_genome_size / 1e9, 2)
            logger.debug('Total genome size set to: {} bp'.format(total_genome_size))
            logger.debug('... human readable: {} Gbp'.format(gs_gbp))
            summary_stats = prepare_summary_statistics(len_stats, total_genome_size, read_count)
            logger.debug('Summary statistics computed')
            os.makedirs(os.path.abspath(os.path.dirname(cargs.summary_output)), exist_ok=True)
            with open(cargs.summary_output, 'w') as stats:
                _ = stats.write('\n'.join(['\t'.join([k, str(v)]) for k, v in summary_stats]))
                _ = stats.write('\n')
            logger.debug('Summary statistics written to file {}'.format(cargs.summary_output))

    logger.debug('Done - processed {} reads'.format(read_count))
    gc_info_reads = sum(gc_stats.values())
    assert gc_info_reads == read_count, \
        'Missing read information: read {} / GC stats {}'.format(read_count, gc_info_reads)

    os.makedirs(os.path.dirname(os.path.abspath(cargs.output)), exist_ok=True)
    dump = {'num_reads': read_count,
            'genome_size': total_genome_size,
            'genome_size_file': cargs.size_file,
            'gc_bins': gc_stats,
            'len_stats': len_stats,
            'nuc_stats': nuc_stats,
            'summary': summary_stats,
            'timestamp': str(ti.ctime())}
    with open(cargs.output, 'wb') as stats_dump:
        _ = pck.dump(dump, stats_dump)
    logger.debug('Statistics saved to output file')

    return


def copy_output_from_existing_source(file_paths, output, logger):
    """
    :param file_paths:
    :param output:
    :param logger:
    :return:
    """
    dest_path = None
    if file_paths is None:
        pass
    else:
        dest_path = None
        for file_path in file_paths:
            if os.path.isfile(file_path):
                # this check is important because it is not guaranteed
                # that other stats computation have finished already
                # (it is likely, though)
                logger.debug('Found existing statistics source at path: {}'.format(file_path))
                os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
                dest_path = sh.copy(file_path, output)
                logger.debug('Copied statistics from {} to {}'.format(file_path, dest_path))
                break
    return dest_path


def main(logger, cargs):
    """
    :param logger:
    :param cargs:
    :return:
    """
    logger.debug('Starting computations')

    logger.debug('Checking if output can be copied from other source...')
    stats_dump_path = copy_output_from_existing_source(cargs.copy_dump, cargs.output, logger)

    if cargs.summary_output is not None:
        _ = copy_output_from_existing_source(cargs.copy_summary, cargs.summary_output, logger)

    if stats_dump_path is None:
        logger.debug('No other sources specified, computing dataset statistics...')
        compute_dataset_statistics(cargs, logger)

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
