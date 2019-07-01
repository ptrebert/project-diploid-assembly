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


FASTQ_QUAL_ENCODING = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""

Read = col.namedtuple('Read', 'line seqid sequence separator qualities record')


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
        help="Full path FASTQ read files (assumed to be gzip or bzip compressed) or BAM alignments.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        dest="output",
        help="Full path to output file (Python pickle format). Directories will"
             " be created if they do not exist.",
    )
    parser.add_argument(
        "--chunk-size",
        '-cs',
        type=int,
        default=5000,
        dest='chunksize',
        help="Read this many reads that are then batch-processed by the worker pool."
             " Default: 5000"
    )
    parser.add_argument(
        "--validate",
        "-val",
        action="store_true",
        default=False,
        dest="validate",
        help="Validate each individual read record in the FASTQ."
             " This increases the runtime substantially. Default: False"
    )
    parser.add_argument(
        "--read-alphabet",
        "-ra",
        type=str,
        default=['A', 'C', 'G', 'T'],
        dest="alphabet",
        help="Specify the alphabet to use for validating the read sequences."
             " Default: A C G T"
    )
    parser.add_argument(
        "--num-cpu",
        "-n",
        type=int,
        default=1,
        dest="numcpu",
        help="Specify number of CPU cores to use for parallel computation. Note"
        " that there is no sanity checking for actually free/available cores."
        " Default: 1",
    )
    return parser.parse_args()


def compute_read_statistics(read):
    """
    :param read:
    :return:
    """
    length = len(read)
    bases = col.Counter(read)
    pct_gc = round((bases['G'] + bases['C']) / length, 3)
    return length, bases, pct_gc


def validate_read_record(check_sequence, check_qualities, read):
    """
    :param read:
    :param check_sequence:
    :param check_qualities:
    :return:
    """
    try:
        assert len(read.seqid) > 0 and read.seqid[0] == '@', \
            'LN {} --- malformed SEQID: {}'.format(read.line, read.seqid)
        assert check_sequence(read.sequence) is not None, \
            'LN {} --- sequence of read {} does not conform to specified alphabet'.format(read.line, read.seqid)
        assert len(read.separator) > 0 and read.separator[0] == '+', \
            'LN {} --- malformed separator line for SEQID: {}'.format(read.line, read.seqid)
        assert check_qualities(read.qualities) is not None, \
            'LN {} --- quality values of read {} do not adhere to FASTQ specification'.format(read.line, read.seqid)
        assert len(read.sequence) == len(read.qualities), \
            'LN {} --- sequence and qualities of different length'.format(read.line)
    except AssertionError as ae:
        ae.args += ('error_record', read.record)
        raise
    return compute_read_statistics(read.sequence)


def read_sequence_records(fpath, chunk_size):
    """
    :param fpath:
    :param chunk_size:
    :return:
    """
    chunk_buffer = [''] * chunk_size
    line_num = 1
    record = 0
    with gz.open(fpath, 'rt', encoding='ascii') as fastq:
        for line in fastq:
            if not line.strip():
                # simply skip over empty lines, though
                # this could be an invalid record
                continue
            if line_num % 4 == 3:
                # every third line ignoring empty lines
                chunk_buffer[record] = line.strip()
                record += 1
                line_num += 1
            else:
                line_num += 1
            if record == chunk_size:
                yield chunk_buffer
                record = 0
                chunk_buffer = [''] * chunk_size
    if record > 0:
        yield chunk_buffer[:record]
    return


def read_complete_records(fpath, chunk_size):
    """
    :param fpath:
    :param chunk_size:
    :return:
    """
    chunk_buffer = [None] * chunk_size
    line_num = 1
    record = 0
    record_is_active = False
    active_record = None
    read_buffer = ['', '', '', '', '', '']
    with gz.open(fpath, 'rt', encoding='ascii') as fastq:
        for ln, line in enumerate(fastq, start=1):
            assert line.endswith('\n'), 'Line {} does not end with newline character'.format(ln)
            if not line.strip() and record_is_active:
                # empty line within active read record
                # is not part of the FASTQ format specification
                raise AssertionError('Empty line (line number {}) within read record {}'.format(ln, active_record))
            elif not line.strip():
                continue
            else:
                if line_num % 4 == 1:
                    assert not record_is_active, \
                        'Encountered new read at line {}, but previous one (SEQID {}) is still active'.format(ln, active_record)
                    record_is_active = True
                    active_record = line.strip()
                    read_buffer[0] = str(ln)
                    read_buffer[1] = active_record
                    line_num += 1
                elif 1 < line_num % 4 < 4:
                    assert record_is_active, \
                        'Inactive record (active SEQID {}) while still reading read information at line {}'.format(active_record, ln)
                    remainder = line_num % 4
                    read_buffer[remainder] = line.strip()
                    line_num += 1
                elif line_num % 4 == 0:
                    assert record_is_active, \
                        'Inactive record (active SEQID {}) while still reading read information at line {}'.format(active_record, ln)
                    read_buffer[4] = line.strip()
                    read_buffer[5] = str(record)
                    chunk_buffer[record] = Read(*read_buffer)
                    record += 1
                    record_is_active = False
                    read_buffer = ['', '', '', '', '', '']
                    line_num += 1
                    if record == chunk_size:
                        yield chunk_buffer
                        record = 0
                        chunk_buffer = [None] * chunk_size
                else:
                    raise RuntimeError('This cannot happen (at line {}: {})'.format(ln, line.strip()))
    if record > 0:
        yield chunk_buffer[:record]
    return


def read_bam_sequence_records(fpath, chunk_size):
    """
    :param fpath:
    :param chunk_size:
    :return:
    """
    chunk_buffer = [''] * chunk_size
    record = 0
    with pysam.AlignmentFile(fpath, 'rb', check_sq=False) as bam:
        for alignment in bam:
            chunk_buffer[record] = alignment.query_sequence
            record += 1
            if record == chunk_size:
                yield chunk_buffer
                chunk_buffer = [''] * chunk_size
                record = 0
    if record > 0:
        yield chunk_buffer[:record]
    return


def dump_error_records(active_chunk, record_num, dump_file):
    """
    :param active_chunk:
    :param record_num:
    :param dump_file:
    :return:
    """
    left_bound = max(0, record_num - 5)
    right_bound = min(len(active_chunk), record_num + 5)
    records_to_write = active_chunk[left_bound:right_bound]
    with gz.open(dump_file, 'wt') as dump:
        for rec in records_to_write:
            _ = dump.write(rec.seqid + '\n')
            _ = dump.write(rec.sequence + '\n')
            _ = dump.write(rec.separator + '\n')
            _ = dump.write(rec.qualities + '\n')
    return


def assemble_file_processors(input_files, validate, logger):
    """
    Quick and dirty - if several input file types should now
    be supported, need to implement an OO approach...

    :param fastq_input:
    :param bam_input:
    :param validate:
    :return:
    """
    seq_re = re.compile("^[" + ''.join(cargs.alphabet) + "]+$")
    qual_re = re.compile("^[" + re.escape(FASTQ_QUAL_ENCODING) + "]+$")
    check_seq = seq_re.match
    check_qual = qual_re.match

    file_paths, chunk_readers, read_processors = [], [], []
    for input_file in input_files:
        file_paths.append(input_file)
        if '.fastq' in input_file and validate:
            chunk_readers.append(read_complete_records)
            read_processors.append(fnt.partial(validate_read_record, *(check_seq, check_qual)))
        elif '.fastq' in input_file:
            chunk_readers.append(read_sequence_records)
            read_processors.append(compute_read_statistics)
        elif '.bam' in input_file:
            logger.debug('Adding BAM processors for file {}'.format(os.path.basename(input_file)))
            chunk_readers.append(read_bam_sequence_records)
            read_processors.append(compute_read_statistics)
        else:
            raise ValueError('Unrecognized input file format: {}'.format(input_file))
    return file_paths, chunk_readers, read_processors


def main(logger, cargs):
    """
    :param logger:
    :param cargs:
    :return:
    """
    logger.debug("Starting computations")

    file_paths, chunk_readers, read_processors = assemble_file_processors(cargs.input, cargs.validate, logger)

    len_stats = col.Counter()
    nuc_stats = col.Counter()
    gc_stats = col.Counter()
    read_count = 0

    gc_bins = np.arange(0, 1.01, 0.01)
    # since numpy.digitize is not called as right-inclusive
    # adjust boundary to include 1
    gc_bins[-1] += 0.001

    gc_buffer = np.zeros(cargs.chunksize, dtype=np.float16)
    with mp.Pool(cargs.numcpu) as pool:
        logger.debug('Initialized worker pool')
        for input_file, chunk_reader, read_processor in zip(file_paths, chunk_readers, read_processors):
            logger.debug('Processing file {}'.format(input_file))
            for chunk in chunk_reader(input_file, cargs.chunksize):
                read_chunk_size = len(chunk)
                logger.debug('Processing chunk of size {}'.format(read_chunk_size))
                try:
                    resit = pool.imap(read_processor, chunk)
                    for idx, (l, n, c) in enumerate(resit):
                        len_stats[l] += 1
                        nuc_stats.update(n)
                        gc_buffer[idx] = c
                    logger.debug('Binning GC values')
                    gc_stats.update(col.Counter(np.digitize(gc_buffer[:read_chunk_size], gc_bins, right=False)))
                    read_count += (len(chunk))
                except AssertionError as ae:
                    if len(ae.args) != 3:
                        raise
                    else:
                        msg = ae.args[0]
                        const = ae.args[1]
                        record_num = int(ae.args[2])
                        logger.error('Found faulty read record: {} - {} - {}'.format(msg, const, record_num))
                        dump_file = os.path.basename(input_file) + 'error'
                        dump_path = os.path.join(os.path.dirname(input_file), dump_file)
                        dump_error_records(chunk, record_num, dump_path)
                        raise

    logger.debug('Done - processed {} reads'.format(read_count))
    gc_info_reads = sum(gc_stats.values())
    assert gc_info_reads == read_count, \
        'Missing read information: read {} / GC stats {}'.format(read_count, gc_info_reads)

    os.makedirs(os.path.dirname(os.path.abspath(cargs.output)), exist_ok=True)
    dump = {'num_reads': read_count,
            'gc_bins': gc_stats,
            'len_stats': len_stats,
            'nuc_stats': nuc_stats,
            'timestamp': str(ti.ctime())}
    with open(cargs.output, 'wb') as stats_dump:
        _ = pck.dump(dump, stats_dump)
    logger.debug('Statistics saved to output file')
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
