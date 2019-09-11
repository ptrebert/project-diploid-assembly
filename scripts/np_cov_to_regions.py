#!/usr/bin/env python3

import os as os
import sys as sys
import logging as log
import io as io
import traceback as trb
import argparse as argp
import collections as col

import numpy as np


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog="np_cov_to_regions.py", description=__doc__)
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
        "--seq-info",
        "-seq",
        type=str,
        required=True,
        dest="seqinfo",
        help="Single line from a FASTA index file (fai).",
    )
    parser.add_argument(
        "--num-regions",
        "-nr",
        type=int,
        required=True,
        dest="numregions",
        help="Number of regions to create",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="stdout",
        help="Output regions to file (or stdout by default)."
    )

    return parser.parse_args()


def read_reference_index(fpath, logger):
    """
    :param fpath:
    :param logger:
    :return:
    """
    logger.debug('Reading sequence info from path {}'.format(fpath))

    with open(fpath, 'r') as fai:
        seq_name, seq_length = fai.readline().split('\t')[:2]

    logger.debug('Read sequence {} (length {}) info file'.format(seq_name, seq_length))

    return seq_name, int(seq_length)


def prepare_data_structures(seqs, logger):
    """
    DEPRECATED

    :param seqs:
    :param logger:
    :return:
    """
    seqs = sorted(seqs, key=lambda x: x[1], reverse=True)
    # debug example: seqs = [('s1', 10), ('s2', 8)]
    total_seq_length = sum([x[1] for x in seqs])
    logger.debug('Maximal total sequence length: {} Gbp'.format(round(total_seq_length / 10**9, 2)))

    pos_depth = np.zeros(total_seq_length, dtype=np.int16)
    seq_boundaries = col.defaultdict(dict)

    last_end = 0
    for sn, sl in seqs:
        offset = last_end - 1
        seq_boundaries[sn]['offset'] = offset
        seq_boundaries[sn]['length'] = sl
        seq_boundaries[sn]['array_start'] = last_end
        seq_boundaries[sn]['array_end'] = last_end + sl
        last_end = last_end + sl

    return seq_boundaries, pos_depth


def read_coverage_per_position(seq_name, seq_length, logger):
    """
    :param seq_name:
    :param seq_length:
    :param logger:
    :return:
    """
    logger.debug('Processing coverage data for sequence {}'.format(seq_name))
    pos_depth = np.zeros(seq_length, dtype=np.int64)

    int64 = np.int64

    for line in sys.stdin:
        _, pos, depth = line.split('\t')
        # -1: bedtools genomecov returns 1-based coords
        pos_depth[int64(pos) - 1] = int64(depth)
    logger.debug('Input stream ended...')
    return pos_depth


def dump_unicov_regions(depth_per_region, pos_depth, seq_name, seq_len, output, logger):
    """
    :param depth_per_region:
    :param pos_depth:
    :param seq_name:
    :param seq_len:
    :param output:
    :param logger:
    :return:
    """
    logger.debug('Writing individual regions to file: {}'.format(output))
    if 'stdout' not in output:
        os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    else:
        output = '/dev/stdout'

    pos_depth = pos_depth.cumsum()

    with open(output, 'w') as dump:
        start = 1
        while 1:
            pos_depth -= depth_per_region
            try:
                current_end = np.nonzero(pos_depth > 0)[0][0]
            except IndexError:
                _ = dump.write('{}:{}-{}\n'.format(seq_name, start, seq_len))
                break
            else:
                _ = dump.write('{}:{}-{}\n'.format(seq_name, start, current_end))
                start = current_end

    logger.debug('Regions with approx. uniform coverage generated')
    return


def main(logger, cargs):
    """
    :param logger:
    :param cargs:
    :return:
    """
    logger.debug("Starting computations")

    seq_name, seq_length = read_reference_index(cargs.seqinfo, logger)
    pos_depth = read_coverage_per_position(seq_name, seq_length, logger)

    total_depth = pos_depth.sum()
    logger.debug('Total cumulative depth: {}'.format(total_depth))

    depth_per_region = np.int64(round(total_depth / cargs.numregions, 0))
    logger.debug('Approx. depth per region ({}): {}'.format(cargs.numregions, depth_per_region))

    dump_unicov_regions(depth_per_region, pos_depth, seq_name, seq_length, cargs.output, logger)
    return


if __name__ == "__main__":
    logger = None
    rc = 0
    try:
        log_msg_format = "%(asctime)s | %(levelname)s | %(message)s"
        cargs = parse_command_line()
        if cargs.debug:
            log.basicConfig(stream=sys.stderr, level=log.DEBUG, format=log_msg_format)
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
