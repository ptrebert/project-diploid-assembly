#!/usr/bin/env python3

import os as os
import sys as sys
import logging as log
import io as io
import traceback as trb
import argparse as argp


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
        "--input-vcf",
        "-iv",
        type=str,
        required=True,
        dest="input",
        help="Full path to input VCF with FLOAT GQ values.",
    )
    parser.add_argument(
        "--output-vcf",
        "-ov",
        type=str,
        required=True,
        dest="output",
        help="Full path to output VCF with INTEGER GQ values.",
    )

    return parser.parse_args()


def read_input_vcf(fpath, logger):
    """
    :param fpath:
    :param logger:
    :return:
    """
    logger.debug('Reading VCF from path {}'.format(fpath))

    out_buffer = io.StringIO()
    with open(fpath, 'r') as vcf_in:
        for line in vcf_in:
            if line.startswith('##FORMAT=<ID=GQ,Number=1,Type=Float,'):
                out_buffer.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
            elif not line.startswith('#'):
                cols = line.strip().split('\t')
                genotype_format = cols[-2].split(':')
                assert genotype_format[1] == 'GQ', 'Unexpected genotype format: {}'.format(line.strip())
                genotype_values = cols[-1].split(':')
                genotype_values[1] = str(int(round(float(genotype_values[1]), 0)))
                genotype_values = ':'.join(genotype_values)
                new_line = '\t'.join(cols[:-1]) + '\t' + genotype_values
                out_buffer.write(new_line + '\n')
            else:
                out_buffer.write(line)

    return out_buffer


def main(logger, cargs):
    """
    :param logger:
    :param cargs:
    :return:
    """
    logger.debug("Starting computations")

    out_buffer = read_input_vcf(cargs.input, logger)

    os.makedirs(os.path.abspath(os.path.dirname(cargs.input)), exist_ok=True)

    logger.debug('Writing normalized VCF to path {}'.format(cargs.output))

    with open(cargs.output, 'w') as vcf_out:
        _ = vcf_out.write(out_buffer.getvalue())

    logger.debug('New VCF dumped')

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
