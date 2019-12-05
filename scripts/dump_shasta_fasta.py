#!/usr/bin/env python3

import os as os
import sys as sys
import logging as log
import io as io
import traceback as trb
import argparse as argp


import dnaio as dnaio


__doc__ = """
Prepare an input FASTA file for the Shasta assembler.
The FASTA file is uncompressed and is not standard-compliant
in the sense that each read (= sequence) is put on a single
line, i.e., no 80 or 120 character limit per line is enforced.
Currently, only FASTQ is supported as input format.
"""


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog="dump_shasta_fasta.py", description=__doc__)
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
        "--input-fq",
        "-i",
        type=str,
        dest="input",
        help="Full path to input FASTQ file (can be compressed)",
    )
    parser.add_argument(
        "--output-fa",
        "-o",
        type=str,
        required=True,
        dest="output",
        help="Full path to output FASTA file. This will NOT be compressed, "
             "even if you specify a file ending like '.gz'",
    )
    parser.add_argument(
        "--buffer-size",
        "-bs",
        type=int,
        default=1,
        dest="buffer_size",
        help="Specify how much memory (approximately) can be used to buffer "
             "reads before dumping them to file. This number is always assumed "
             "to be gigabyte. Default: 1 (gigabyte)"
    )

    return parser.parse_args()


def main(logger, cargs):
    """
    :param logger:
    :param cargs:
    :return:
    """
    logger.debug("Start processing FASTQ input from: {}".format(cargs.input))
    if not os.path.isfile(cargs.input):
        logger.error('The input path you specified seems not to be valid: {}'.format(cargs.input))
        raise FileNotFoundError(cargs.input)

    if cargs.buffer_size < 1:
        raise ValueError('Buffer size (in GB) has to be a value >= 1')

    input_info = os.stat(os.path.abspath(cargs.input))
    input_size = input_info.st_size
    logger.debug('Input file size: {} B'.format(input_size))
    logger.debug('Input file size: {} GB'.format(round(input_size / (1024 ** 3), 2)))

    buffer_limit = cargs.buffer_size * (1000 ** 3)  # GB -> MB -> KB -> B

    out_buffer = io.StringIO()
    buffer_size = sys.getsizeof('')
    num_records = 0
    first_dump = True

    total_written = 0

    with dnaio.open(cargs.input) as fastq:
        for read in fastq:
            num_records += 1
            buffer_size += len(read) + len(read.name)
            out_buffer.write('>' + read.name + '\n' + read.sequence + '\n')

            if buffer_size >= buffer_limit:
                open_mode = 'w' if first_dump else 'a'
                if first_dump:
                    logger.debug('Creating FASTA output path at: {}'.format(cargs.output))
                    os.makedirs(os.path.dirname(os.path.abspath(cargs.output)), exist_ok=True)
                    logger.debug('Path created')
                first_dump = False

                with open(cargs.output, open_mode) as fasta:
                    _ = fasta.write(out_buffer.getvalue())

                logger.debug('Buffer dumped to output')
                out_buffer = io.StringIO()
                total_written += buffer_size
                logger.debug('Records processed: {}'.format(num_records))
                logger.debug('Approx. bytes dumped to output: {} B'.format(total_written))
                buffer_size = sys.getsizeof('')

    open_mode = 'w' if first_dump else 'a'
    if first_dump:
        logger.debug('Creating FASTA output path at: {}'.format(cargs.output))
        os.makedirs(os.path.dirname(os.path.abspath(cargs.output)), exist_ok=True)
        logger.debug('Path created')

    with open(cargs.output, open_mode) as fasta:
        _ = fasta.write(out_buffer.getvalue())

    logger.debug('Buffer dumped to output')
    total_written += buffer_size

    logger.debug('Total records processed: {}'.format(num_records))
    logger.debug('Approx. bytes dumped to output: {} B'.format(total_written))
    logger.debug('Approx. GB dumped to output: {} GB'.format(round(total_written / (1023 ** 3), 2)))

    output_info = os.stat(os.path.abspath(cargs.output))
    output_size = output_info.st_size
    logger.debug('Output file size: {} B'.format(output_size))
    logger.debug('Output file size: {} GB'.format(round(output_size / (1024 ** 3), 2)))

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
