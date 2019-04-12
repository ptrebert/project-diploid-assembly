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
import collections as col
import multiprocessing as mp


import numpy as np


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog="get_stranded_read_counts.py", description=__doc__)
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
        "--fastq-input",
        "-fi",
        type=str,
        required=True,
        dest="fastq",
        nargs='+',
        help="Full path FASTQ read files (assumed to be gzip or bzip compressed).",
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


def main(logger, cargs):
    """
    :param logger:
    :param cargs:
    :return:
    """
    logger.debug("Starting computations")

    len_stats = col.Counter()
    nuc_stats = col.Counter()
    gc_stats = col.Counter()
    read_count = 0

    gc_bins = np.arange(0, 1.01, 0.01)
    # since numpy.digitize is not called as right-inclusive
    # adjust boundary to include 1
    gc_bins[-1] += 0.001

    read_buffer = [''] * cargs.chunksize
    gc_buffer = np.zeros(cargs.chunksize, dtype=np.float16)
    r = 0
    with mp.Pool(cargs.numcpu) as pool:
        logger.debug('Initialized worker pool')
        for fastq in cargs.fastq:
            with gz.open(fastq, mode='rt') as readfile:
                for line in readfile:
                    if line[0] in ['A', 'C', 'G', 'T']:
                        read_buffer[r] = line.strip()
                        r += 1
                    if r == cargs.chunksize:
                        resit = pool.imap_unordered(compute_read_statistics, read_buffer)
                        logger.debug('Processing chunk')
                        for idx, (l, n, c) in enumerate(resit):
                            len_stats[l] += 1
                            nuc_stats.update(n)
                            gc_buffer[idx] = c
                        logger.debug('Binning GC values')
                        gc_stats.update(col.Counter(np.digitize(gc_buffer, gc_bins, right=False)))
                        logger.debug('Resetting read counter')
                        read_count += r
                        r = 0

        logger.debug('Processing last chunk of size {}'.format(r))
        read_buffer = read_buffer[:r]
        gc_buffer = gc_buffer[:r]
        read_count += r
        resit = pool.imap_unordered(compute_read_statistics, read_buffer)
        logger.debug('Processing chunk')
        for idx, (l, n, c) in enumerate(resit):
            len_stats[l] += 1
            nuc_stats.update(n)
            gc_buffer[idx] = c
        logger.debug('Binning GC values')
        gc_stats.update(col.Counter(np.digitize(gc_buffer, gc_bins, right=False)))

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
