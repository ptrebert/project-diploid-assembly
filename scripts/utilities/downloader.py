#!/usr/bin/env python3

import os as os
import sys as sys
import stat as stat
import logging as log
import io as io
import traceback as trb
import argparse as argp
import subprocess as sp

CMD_DL_COMPRESSED_PARALLEL = 'aria2c --out={output} --file-allocation=none -s {num_conn} -x {num_conn} {remote_path}'

CMD_DL_COMPRESSED_SINGLE = 'wget --no-verbose -O {output} {remote_path}'

CMD_DL_UNCOMPRESSED_SINGLE = 'wget --no-verbose -O /dev/stdout {remote_path} | gzip > {output}'

CMD_COPY_LOCAL = 'cp {remote_path} {output}'

CMD_SYM_LINK = 'ln -s {remote_path} {output}'


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

    mut_group = parser.add_mutually_exclusive_group()
    mut_group.add_argument(
        "--request-file",
        "-req",
        type=str,
        default=None,
        dest="request_file",
        help="Process generic request file."
    )
    mut_group.add_argument(
        "--ena-file-report",
        "-ena",
        type=str,
        default=None,
        dest="ena_file_report",
        help="Download ENA file report/metadata."
    )
    mut_group.add_argument(
        "--shasta-exec",
        "-sxc",
        type=str,
        default=None,
        dest="shasta_exec",
        help="Download Shasta executable."
    )

    parser.add_argument(
        "--parallel-conn",
        "-p",
        type=int,
        default=1,
        dest="parallel_conn",
        help="If the download happens via aria2c, allow this many connections in parallel."
    )

    parser.add_argument(
        "--force-local-copy",
        "-flc",
        action="store_true",
        default=False,
        dest="force_local_copy",
        help="For files located on localhost, make a copy instead of soft linking. Default: False"
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        dest="output",
        help="Path to output file (local destination of download)."
    )

    shasta_grp = parser.add_argument_group('Shasta download arguments')

    shasta_grp.add_argument(
        "--shasta-version",
        "-sver",
        type=str,
        dest="shasta_version",
        help="Shasta version to check."
    )

    shasta_grp.add_argument(
        "--shasta-path",
        "-spth",
        type=str,
        dest="shasta_path",
        help="Local path for Shasta executable. This is assumed to be '$CONDA_PREFIX/bin/shasta'."
    )

    return parser.parse_args()


def handle_request_file_download(req_file_path, output_path, parallel_conn, force_copy, logger):
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
        local_path = req_file.readline().strip()

    if local_path != output_path:
        logger.error('Path mismatch between output and request file: {} / {}'.format(output_path, local_path))
        raise ValueError('Path mismatch between output and request file: {} / {}'.format(output_path, local_path))

    logger.debug('Creating folder hierarchy')
    local_folder = os.path.dirname(os.path.abspath(local_path))
    os.makedirs(local_folder, exist_ok=True)

    if os.path.isfile(remote_path):
        logger.debug('Request file contains local resource')
        if force_copy:
            selected_method = CMD_COPY_LOCAL
        else:
            common_prefix = os.path.commonprefix([remote_path, local_folder])
            if common_prefix == '/':
                logger.warning('Creating sym link between different file systems')
            selected_method = CMD_SYM_LINK

        call = selected_method.format(
                **{
                    'remote_path': remote_path,
                    'output': output_path
                })

    elif any([remote_path.endswith(x) for x in ['.gz', '.bam', '.bam.1']]):
        logger.debug('Downloading remote path in parallel with {} connections'.format(parallel_conn))
        call = CMD_DL_COMPRESSED_PARALLEL.format(
            **{
                'remote_path': remote_path,
                'num_conn': parallel_conn,
                'output': output_path
            })
    else:
        logger.debug('Downloading remote path with single process')
        call = CMD_DL_UNCOMPRESSED_SINGLE.format(
            **{
                'remote_path': remote_path,
                'output': output_path
            })

    logger.debug('Executing system call: {}'.format(call))

    try:
        out = sp.check_output(call, shell=True, env=None, stderr=sp.STDOUT)
        logger.debug('System call completed... dumping output to logfile...')
        logger.debug('=== begin: stdout/err')
        logger.debug(out.decode('utf-8'))
        logger.debug('=== end: stdout/err')
    except sp.CalledProcessError as spe:
        logger.error('System call failed with exit code: {}'.format(spe.returncode))
        logger.error('Error output: ---')
        logger.error(spe.output)
        raise spe

    return


def handle_ena_file_report_download(ena_dl_url, output_path, logger):
    """
    :param ena_dl_url:
    :param output_path:
    :param logger:
    :return:
    """
    logger.debug('Creating folder hierarchy')
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    call = CMD_DL_COMPRESSED_SINGLE.format(
        **{
            'remote_path': ena_dl_url,
            'output': output_path
        }
    )
    logger.debug('Executing system call: {}'.format(call))

    try:
        out = sp.check_output(call, shell=True, env=None, stderr=sp.STDOUT)
        logger.debug('System call completed... dumping output to logfile...')
        logger.debug('=== begin: stdout/err')
        logger.debug(out.decode('utf-8'))
        logger.debug('=== end: stdout/err')
    except sp.CalledProcessError as spe:
        logger.error('System call failed with exit code: {}'.format(spe.returncode))
        logger.error('Error output: ---')
        logger.error(spe.output)
        raise spe

    return


def handle_shasta_executable_download(shasta_dl_url, shasta_local_path, shasta_version, shasta_ver_check, logger):
    """
    :param shasta_dl_url:
    :param shasta_local_path:
    :param shasta_version:
    :param shasta_ver_check:
    :param logger:
    :return:
    """

    conda_path = None
    for key, value in dict(os.environ).items():
        if 'conda' in key.lower():
            if key == 'CONDA_PREFIX':
                conda_path = value
            else:
                logger.debug('Skipping over conda ENV VAR: {} - {}'.format(key, value))
        if 'path' in key.lower():
            logger.debug('INFO PATH - {}: {}'.format(key, value))
    if conda_path is None:
        logger.write('Error: could not identify conda environment under '
                     'environment variable CONDA_PREFIX. Shasta may not '
                     'end up in the desired Conda environment. Aborting...')
        raise ValueError('No conda path detected')

    conda_path = os.path.join(conda_path, 'bin')

    if conda_path not in shasta_local_path:
        logger.error('Detected accessible Conda path is different from specified '
                     'output Conda path - aborting: {} vs {}'.format(conda_path, shasta_local_path))
        raise ValueError('Conda path mismatch: {} vs {}'.format(conda_path, shasta_local_path))

    if shasta_local_path.strip('/').endswith('/bin'):
        logger.warning('I dare to modify the specified output path because I think '
                       'the name of the Shasta executable is missing.\nOriginal: {}\n'
                       'Modified: {}'.format(shasta_local_path, os.path.join(shasta_local_path, 'shasta')))
        shasta_local_path = os.path.join(shasta_local_path, 'shasta')

    if not shasta_local_path.endswith('shasta'):
        logger.error('Specified local path for Shasta executable does not end in "shasta": {}'.format(shasta_local_path))
        raise ValueError('Local Shasta path does not end in "shasta".')

    if not os.path.isdir(os.path.dirname(shasta_local_path)):
        logger.error('The specified output path (should be "/bin" in Conda environment) does not exist. '
                     'I will not create that path because the Conda environment '
                     'should already exist: {}'.format(os.path.dirname(shasta_local_path)))
        raise RuntimeError('Conda environment seems not to exist: {}'.format(os.path.dirname(shasta_local_path)))

    download_executable = True
    if os.path.isfile(shasta_local_path):
        logger.debug('Existing shasta executable found - removing...')
        try:
            os.unlink(shasta_local_path)
        except (OSError, IOError) as error:
            logger.warning('WARNING: could not remove existing shasta executable at path {}'.format(shasta_local_path))
            logger.error('ERROR message: {}'.format(str(error)))
            logger.warning('Proceeding - will fail if existing shasta version does not match...')
            download_executable = False
    if download_executable:
        logger.debug('Placing shasta executable at path: {}'.format(shasta_local_path))
        call = CMD_DL_COMPRESSED_SINGLE.format(
            **{
                'remote_path': shasta_dl_url,
                'output': shasta_local_path
            }
        )
        logger.debug('Executing system call: {}'.format(call))

        try:
            out = sp.check_output(call, shell=True, env=None, stderr=sp.STDOUT)
            logger.debug('System call completed... dumping output to logfile...')
            logger.debug('=== begin: stdout/err')
            logger.debug(out.decode('utf-8'))
            logger.debug('=== end: stdout/err')
        except sp.CalledProcessError as spe:
            logger.error('System call failed with exit code: {}'.format(spe.returncode))
            logger.error('Error output: ---')
            logger.error(spe.output)
            raise spe
    logger.debug('Attempt of changing permissions to user-rwx for Shasta executable...')
    os.chmod(shasta_local_path, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
    logger.debug('Success')
    try:
        shasta_ver = sp.check_output('shasta --version',
                                     stderr=sp.STDOUT,
                                     shell=True)
        shasta_version_info = shasta_ver.decode('utf-8')
        if shasta_version in shasta_version_info:
            logger.debug('Shasta assembler available in PATH, version confirmed...')
            with open(shasta_ver_check, 'w') as dump:
                _ = dump.write(shasta_version_info)
        else:
            logger.error('Shasta version mismatch - found info: '
                         '{}\nRequired: {}'.format(shasta_version_info, shasta_version))
            raise ValueError('Incompatible Shasta version '
                             '({} required): {}'.format(shasta_version, shasta_version_info))
    except sp.CalledProcessError as spe:
        rc = spe.returncode
        err_msg = str(spe.output)
        logger.error('Shasta version check call returned {}: {}'.format(rc, err_msg))
        raise spe
    return


def main(logger, cargs):
    """
    :param logger:
    :param cargs:
    :return:
    """
    if cargs.request_file is not None:
        logger.debug('Handling request file download: {}'.format(cargs.request_file))
        handle_request_file_download(
            cargs.request_file,
            cargs.output,
            cargs.parallel_conn,
            cargs.force_local_copy,
            logger
        )
    elif cargs.ena_file_report is not None:
        logger.debug('Handling ENA file report/metadata request')
        handle_ena_file_report_download(
            cargs.ena_file_report,
            cargs.output,
            logger
        )
    elif cargs.shasta_exec is not None:
        logger.debug('Downloading Shasta executable')
        handle_shasta_executable_download(
            cargs.shasta_exec,
            cargs.shasta_path,
            cargs.shasta_version,
            cargs.output,
            logger
        )
    else:
        raise RuntimeError('No valid download option selected')

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
