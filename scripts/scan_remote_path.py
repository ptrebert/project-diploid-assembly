#!/usr/bin/env python3

import os as os
import sys as sys
import logging as log
import io as io
import json as json
import re as re
import traceback as trb
import argparse as argp
import collections as col
import ftplib as ftplib


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog="scan_remote_path.py", description=__doc__)
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
        "--server",
        "-srv",
        type=str,
        default="ftp.1000genomes.ebi.ac.uk",
        dest="server",
        help="Remote server URL. Default: ftp.1000genomes.ebi.ac.uk",
    )
    parser.add_argument(
        "--ftp-path",
        "-fp",
        type=str,
        required=True,
        dest="ftp_path",
        help="Folder on FTP server to scan for input",
    )
    parser.add_argument(
        "--collect-files",
        "-cf",
        type=str,
        nargs='+',
        dest="collect_files",
        help="Specify file extensions to collect from remote path."
    )
    parser.add_argument(
        "--sort-files",
        "-sf",
        type=str,
        nargs='+',
        dest="sort_files",
        help="For each file extension, specify the local path prefix for the file."
    )
    parser.add_argument(
        "--assume-pacbio-native",
        "-pbn",
        action="store_true",
        default=False,
        dest="pacbio_native",
        help="Assume that BAM files are in PacBio-native format (extension: pbn.bam)"
    )
    parser.add_argument(
        "--file-infix",
        "-in",
        type=str,
        required=True,
        dest="infix",
        help="Define the infix for the output files."
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        dest="output",
        help="Path to output JSON file"
    )

    return parser.parse_args()


def annotate_remote_files(remote_files, cargs, logger):
    """
    :param remote_files:
    :param cargs:
    :param logger:
    :return:
    """
    meta_files = ['readme', 'manifest']
    match_individual = re.compile('[A-Z0-9]+')
    file_collector = col.defaultdict(list)

    for full_path in remote_files:
        path, name = os.path.split(full_path)
        if any([x in name.lower() for x in meta_files]):
            continue
        if not any([name.endswith(x) for x in cargs.collect_files]):
            continue
        mobj = match_individual.match(name)
        if mobj is None:
            logger.warning('No individual identified for file: {}'.format(name))
            continue
        individual = mobj.group(0)
        logger.debug('Extracted individual {} for file {}'.format(individual, name))
        if any([x in name.lower() for x in ['ccs', 'q20']]):
            tech = 'ccs'
        elif any([x in name.lower() for x in ['clr', 'subreads']]):
            tech = 'clr'
        else:
            logger.warning('Skipping file {} - could not id seq. tech.'.format(name))
            continue
        local_path = None
        file_ext = None
        for lp, ext in zip(cargs.sort_files, cargs.collect_files):
            if name.endswith(ext):
                local_path = lp
                file_ext = ext
                break
        if local_path is None:
            raise ValueError('Could not identify extension / local path for file: {}'.format(name))
        if file_ext.lower() == 'bam' and cargs.pacbio_native:
            file_ext = 'pbn.bam'
        elif file_ext.lower() == 'bam':
            file_ext = 'sam.bam'
        else:
            pass
        base_file_prefix = '_'.join([individual, cargs.infix + tech])
        file_collector[(local_path, base_file_prefix, file_ext)].append(full_path)
    logger.debug('Identified {} different file groups'.format(len(file_collector)))
    return file_collector


def enumerate_file_parts(file_groups, cargs, logger):
    """
    :param file_groups:
    :param cargs:
    :param logger:
    :return:
    """
    output = dict()
    for (local_path, file_prefix, file_ext), remote_files in file_groups.items():
        logger.debug('Processing {} / {}'.format(file_prefix, file_ext))
        if len(remote_files) == 1:
            raise ValueError('Single remote file - "complete" suffix needed')
        for part_num, rf in enumerate(remote_files, start=1):
            key = os.path.join(local_path, file_prefix + '.part' + str(part_num))
            full_local_path = key + '.' + file_ext
            assert key not in output, 'Duplicate key: {}\n\n{}\n'.format(key, file_groups)
            remote_path = os.path.join(cargs.server, rf)
            if not remote_path.startswith('ftp://'):
                remote_path = 'ftp://' + remote_path
            output[key] = {
                "remote_path": remote_path,
                "local_path": full_local_path
            }
    logger.debug('Prepared annotation for {} file splits'.format(len(output)))
    return output


def main(logger, cargs):
    """
    :param logger:
    :param cargs:
    :return:
    """
    if not len(cargs.collect_files) == len(cargs.sort_files):
        raise ValueError('Need one path prefix per file extension (sort and collect files)')
    logger.debug("Starting remote FTP scan...")
    server = ftplib.FTP(cargs.server)
    logger.debug('Attempt for anonymous login')
    server.login()

    remote_files = server.nlst(cargs.ftp_path)
    logger.debug("Identified {} remote files under path {}".format(len(remote_files), cargs.ftp_path))

    file_groups = annotate_remote_files(remote_files, cargs, logger)

    output = enumerate_file_parts(file_groups, cargs, logger)

    os.makedirs(os.path.dirname(os.path.abspath(cargs.output)), exist_ok=True)
    with open(cargs.output, 'w') as dump:
        json.dump(output, dump, sort_keys=True, indent=2)

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
