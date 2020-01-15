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


CCS_TECH_KEYWORDS = ['ccs', 'q20']
CLR_TECH_KEYWORDS = ['clr']


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
        "--assume-clr-subreads",
        "-clr",
        action="store_true",
        default=False,
        dest="clr_subreads",
        help="If no CCS indicator ({}) is part of the filename, and the file is a "
             "subreads file, then assume that it is a PacBio CLR dataset".format(CCS_TECH_KEYWORDS)
    )
    parser.add_argument(
        "--assume-paired-reads",
        "-prd",
        action="store_true",
        default=False,
        dest="paired_reads",
        help="Assume paired reads with mate indicators R1/R2"
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
        "--file-suffix",
        "-sfx",
        type=str,
        default=None,
        dest="suffix",
        help="Append suffix to all file names (before part or mate indicator)."
    )
    parser.add_argument(
        "--fix-tech",
        "-ft",
        type=str,
        default=None,
        dest="fix_tech",
        help="Use a fix technology for all files irrespective of tech indicators in the file name."
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        dest="output",
        help="Path to output JSON file"
    )

    return parser.parse_args()


def extract_maximal_match(file_name):
    """
    :param file_name:
    :return:
    """
    lib_id = re.compile('[A-Za-z0-9_\-]+')
    longest_match = 0
    matched_string = ''
    for match in re.finditer(lib_id, file_name):
        current_match = match.group(0)
        if len(current_match) > longest_match:
            longest_match = len(current_match)
            matched_string = current_match
    if longest_match == 0:
        raise ValueError('Could not match name/library identifier in file name: {}'.format(file_name))
    return matched_string


def collect_strong_tech_indicators(remote_files, skip_files):
    """
    In case there is a tech mix in one folder, try to separate
    different techs by checking if there is a library identifier
    (typically a long, convoluted string with date etc.) that is
    tied to a strong tech indicator (e.g., CCS, Q20 etc). If so,
    that library identifier has to be ignored in case of potentially
    ambiguous matches for, e.g., "subreads" datasets.

    :param remote_files:
    :param skip_files:
    :return:
    """
    lib_tech_lut = dict()
    for full_path in remote_files:
        path, name = os.path.split(full_path)
        if any([x in name.lower() for x in skip_files]):
            continue
        library_id = extract_maximal_match(name)
        if any([x in name.lower() for x in CCS_TECH_KEYWORDS]):
            tech = 'ccs'
        elif any([x in name.lower() for x in CLR_TECH_KEYWORDS]):
            tech = 'clr'
        else:
            tech = None
        if library_id in lib_tech_lut:
            known_tech = lib_tech_lut[library_id]
            if known_tech is not None:
                if tech is not None:
                    if known_tech != tech:
                        raise ValueError('Library ID {} has two tech '
                                         'indicators: {} / {}'.format(library_id, known_tech, tech))
                else:
                    # strong tech indicator exists for lib
                    # this is fine
                    pass
            else:
                lib_tech_lut[library_id] = tech
        else:
            lib_tech_lut[library_id] = tech
    return lib_tech_lut


def annotate_remote_files(remote_files, cargs, logger):
    """
    :param remote_files:
    :param cargs:
    :param logger:
    :return:
    """
    meta_files = ['readme', 'manifest']
    aux_files = ['_stats', 'scraps', 'subreadset']
    match_individual = re.compile('[A-Z0-9]+')
    file_collector = col.defaultdict(list)

    logger.debug('Collecting strong tech indicators per library')
    tech_collector = collect_strong_tech_indicators(remote_files, meta_files + aux_files)
    libraries = set(tech_collector.keys())
    logger.debug('Identified {} libraries in remote files'.format(len(tech_collector)))

    for full_path in remote_files:
        path, name = os.path.split(full_path)
        if any([x in name.lower() for x in meta_files]):
            logger.debug('>>> Skipping file (metadata): {}'.format(name))
            continue
        if any([x in name.lower() for x in aux_files]):
            logger.debug('>>> Skipping file (tech./aux.): {}'.format(name))
            continue
        if not any([name.endswith(x) for x in cargs.collect_files]):
            logger.debug('>>> Skipping file b/c of unmatched file extension: {}'.format(name))
            continue
        mobj = match_individual.match(name)
        if mobj is None:
            logger.warning('No individual identified for file: {}'.format(name))
            continue
        individual = mobj.group(0)
        tech = None
        logger.debug('Extracted individual {} for file {}'.format(individual, name))
        if cargs.fix_tech is not None:
            tech = cargs.fix_tech
        elif any([x in name.lower() for x in CCS_TECH_KEYWORDS]):
            tech = 'ccs'
        elif any([x in name.lower() for x in CLR_TECH_KEYWORDS]):
            tech = 'clr'
        elif tech is None and 'subreads' in name.lower() and cargs.clr_subreads:
            # need to check if library ID is associated with a strong tech indicator
            for lib in libraries:
                if lib in name:
                    tech = tech_collector[lib]
                    if tech is None:
                        logger.warning('Assuming CLR/subreads file: {}'.format(name))
                        tech = 'clr'
                    break
        else:
            logger.warning('>>> Skipping file {} - could not id seq. tech.'.format(name))
            continue
        if tech == 'ccs' and 'subreads' in name.lower():
            logger.debug('>>> Skipping file (CCS/subreads): {}'.format(name))
            continue
        logger.debug('Accepted file: {}'.format(name))
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
        if cargs.suffix is None:
            base_file_prefix = '_'.join([individual, cargs.infix + tech])
        else:
            base_file_prefix = '_'.join([individual, cargs.infix + tech, cargs.suffix])
        logger.debug('Adding file to collection: {}'.format(name))
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
            if cargs.paired_reads:
                if part_num > 2:
                    raise ValueError('Assuming paired reads, but mate indicator > 2: {}'.format(remote_files))
                key = os.path.join(local_path, file_prefix + '.R' + str(part_num))
            else:
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
