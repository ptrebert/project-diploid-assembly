#!/usr/bin/env python3

import os as os
import sys as sys
import io as io
import glob as glob
import traceback as trb
import argparse as argp
import subprocess as sp


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog="file_reset.py", description=__doc__)
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
        "--request-file",
        "-req",
        type=str,
        default=None,
        dest="request_file",
        help="Process generic request file."
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        dest="output",
        help="Path to output file (local destination of download)."
    )

    return parser.parse_args()


def handle_request_file(req_file_path, output_path, debug):
    """
    :param req_file_path:
    :param output_path:
    :param parallel_conn:
    :param force_copy:
    :param logger:
    :return:
    """

    with open(req_file_path, 'r') as req_file:
        remote_path = req_file.readline().strip()
        remote_name = os.path.basename(remote_path)
        local_path = req_file.readline().strip()
        local_name = os.path.basename(local_path)
    
    local_path = os.path.split(os.path.dirname(req_file_path))[0]

    input_file_path = os.path.join(local_path, local_name)
    output_file_path = os.path.join(output_path, remote_name)

    call = 'rsync {} {}'.format(input_file_path, output_file_path)
    if debug:
        print('Would execute:\n{}\n'.format(call))
        return

    try:
        out = sp.check_output(call, shell=True, env=None, stderr=sp.STDOUT)
    except sp.CalledProcessError as spe:
        sys.stderr.write('\nSystem call failed with exit code: {}'.format(spe.returncode))
        sys.stderr.write('\nError output: ---')
        sys.stderr.write('\n>>>\n' + spe.output.decode('utf-8').strip() + '\n')
        raise spe

    return


def main(cargs):
    """
    :param cargs:
    :return:
    """
    req_files = []
    if os.path.isfile(cargs.request_file):
        req_files = [cargs.request_file]
    elif os.path.isdir(cargs.request_file):
        req_files = glob.glob(os.path.join(cargs.request_file, '*.request'))
    else:
        raise RuntimeError('Cannot handle path: {}'.format(cargs.request_file))

    os.makedirs(os.path.abspath(cargs.output), exist_ok=True)

    for rf in req_files:
        print('Processing {}'.format(rf))
        handle_request_file(rf, os.path.abspath(cargs.output), cargs.debug)

    return


if __name__ == "__main__":
    logger = None
    rc = 0
    try:
        cargs = parse_command_line()
        main(cargs)
    except Exception as exc:
        rc = 1
        trb.print_exc()
    finally:
        sys.exit(rc)
