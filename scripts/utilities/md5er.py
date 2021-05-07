#!/usr/bin/env python3

import os as os
import sys as sys
import io as io
import glob as glob
import traceback as trb
import argparse as argp
import subprocess as sp
import multiprocessing as mp


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog="md5er.py", description=__doc__)
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
        "--folder",
        "-f",
        type=str,
        dest="folder",
        help="Path to data folder."
    )

    parser.add_argument(
        '--jobs',
        '-j',
        type=int,
        dest='jobs',
        help='Number of parallel jobs'
    )

    parser.add_argument(
        '--file-ext',
        '-x',
        type=str,
        dest='file_ext',
        default='',
        help='Select only input files having this file extension.'
    )

    return parser.parse_args()


def compute_md5(file_path):

    out_path = file_path + '.md5'
    
    call = 'md5sum {} > {}'.format(file_path, out_path)

    try:
        out = sp.check_output(call, shell=True, env=None, stderr=sp.STDOUT)
    except sp.CalledProcessError as spe:
        sys.stderr.write('\nSystem call failed with exit code: {}'.format(spe.returncode))
        sys.stderr.write('\nError output: ---')
        sys.stderr.write('\n>>>\n' + spe.output.decode('utf-8').strip() + '\n')
        raise spe
    return out_path


def main(cargs):
    """
    :param cargs:
    :return:
    """
    all_selected_files = []
    for root, dirs, files in os.walk(cargs.folder, followlinks=False):
        selected_files = [os.path.join(root, f) for f in files if f.endswith(cargs.file_ext)]
        all_selected_files.extend(selected_files)

    with mp.Pool(cargs.jobs) as pool:
        resit = pool.imap_unordered(compute_md5, all_selected_files)
        for res in resit:
            print('Done: ', res)

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
