#!/usr/bin/env python3

import os
import sys
import re
import time
import random
import logging
from logging.handlers import RotatingFileHandler
import subprocess as sp

random.seed()

FILE_PATH_STATUS_LOG = None

# Change here: frequent I/O hiccups on HILBERT
# leading to a FileNotFoundError for abspath.
# Quite annoying...
try:
    script_path = os.path.abspath(__file__)
except FileNotFoundError:
    time.sleep(0.5)  # not sure if Snakemake has a timeout on the status call
    try:
        script_path = os.path.abspath(__file__)
    except FileNotFoundError:
        # to avoid an error condition because of an temp. I/O
        # problem, tell Snakemake everything is fine
        sys.stdout.write('{}\n'.format('running'))
        sys.exit(0)


if FILE_PATH_STATUS_LOG is None:
    FILE_PATH_STATUS_LOG = os.path.join(
        os.path.dirname(
            os.path.abspath(script_path)),
        'log', 'cluster_status.log'
    )

log_msg_format = "%(asctime)s | %(levelname)s | %(message)s"
cluster_status_log = FILE_PATH_STATUS_LOG
os.makedirs(os.path.dirname(cluster_status_log), exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format=log_msg_format,
    handlers=[
        RotatingFileHandler(
            cluster_status_log,
            maxBytes=5000000,
            backupCount=10,
            mode='a'
        )
    ]
)

logger = logging.getLogger()


def parse_qstat_output(job_info, job_id):
    """
    :param job_info:
    :param job_id:
    :return:
    """
    keep_info = [
        'Job_Name',
        'resources_used.mem',
        'resources_used.ncpus',
        'resources_used.walltime',
        'job_state',
        'queue',
        'exec_host',
        'exec_vnode',
        'ctime',
        'mtime',
        'qtime',
        'Resource_List.mem',
        'Resource_List.ncpus',
        'Resource_List.select',
        'Resource_List.walltime',
        'stime',
        'session_id',
        'etime',
        'Exit_status',
        ]

    job_state_codes = {
        'B': 'running',
        'E': 'done',
        'F': 'done',
        'H': 'running',
        'M': 'running',
        'Q': 'running',
        'R': 'running',
        'S': 'running',
        'T': 'running',
        'U': 'running',
        'W': 'running',
        'X': 'done'}

    log_info, job_status, exit_code = '\n=====\nJOB: {}\n'.format(job_id), None, None
    for line in job_info.split('\n'):
        line = line.strip()
        if not line:
            continue
        if 'Job Id:' in line:
            log_info += 'Processing status for job: {}\n'.format(line.split()[-1])
        elif any([x in line for x in keep_info]):
            log_info += '... ' + line + '\n'
            if 'job_state' in line:
                state_code = line.split()[-1]
                try:
                    simple_job_state = job_state_codes[state_code]
                except KeyError as ke:
                    logger.error('Job state lookup failed: {} / {}'.format(state_code, line))
                    raise ke
                else:
                    job_status = simple_job_state
                    log_info += '... Job determined as ' + job_status + '\n'
            if 'Exit_status' in line:
                exit_code = int(line.split()[-1])
                if 0 < exit_code <= 128:
                    job_status = 'failed'
                elif exit_code == 0:
                    job_status = 'success'
                elif exit_code == -3:  # implicit: job status is running
                    # special PBS_VALUE: -3 = JOB_EXEC_RETRY : Job exec failed, do retry
                    # Job could not be executed, but will be retried
                    job_status = 'running'
                elif exit_code < 0:
                    # https://www.pbsworks.com/pdfs/PBSAdminGuide18.2.pdf
                    job_status = 'failed'
                    log_info += 'ERROR: job has negative exit status, indicates that job could not be started\n'
                elif exit_code == 271:
                    log_info += 'ERROR: job has exit code 271 - killed by PBS scheduler for exceeding ' \
                                'CPU or memory resources or time limit.\n'
                    job_status = 'failed'
                elif exit_code > 128:
                    log_info += 'ERROR: job has exit code {} (> 128), ' \
                                'likely killed by signal {}\n'.format(exit_code, exit_code - 128)
                    job_status = 'failed'
                else:
                    log_info += 'ERROR: job has non-zero exit code, and status done: failed for unknown reason\n'
                    job_status = 'failed'

        else:
            continue
    if job_status == 'failed':
        logger.error(log_info)
    elif job_status == 'success':
        logger.info(log_info)
    else:
        if random.randint(0, 100) < 5:
            # for running/ongoing jobs, only sporadically log status
            logger.info(log_info)

    if job_status == 'running' and exit_code is not None and exit_code != -3:
        logger.error('Job {} determined as running, but has error exit code: {}'.format(job_id, exit_code))
        logger.error(log_info)
        job_status = 'failed'

    if job_status == 'done':
        # does not yet have an exit code
        # maybe caught in completing... wait
        job_status = 'running'

    return job_status


job_id = sys.argv[1]

# logger.info('>>> Received input job id: {}'.format(job_id))

mobj = re.match('[0-9]+', job_id)

if mobj is None:
    logger.error('Could not determine numeric part of job ID - exiting...')
    raise ValueError('Malformed job id')

num_job_id = mobj.group(0)
report_job_status = 'failed'
try:
    qstat_output = sp.check_output(
        'qstat -xf {}'.format(num_job_id),
        stderr=sp.STDOUT,
        shell=True,
        timeout=60
    )
except sp.CalledProcessError as cpe:
    if cpe.returncode == 153:
        # job id no longer known to qstat, assume it was fine
        # Snakemake will check output files
        report_job_status = 'success'
    else:
        logger.error('qstat call failed: {} / {}'.format(cpe.returncode, cpe.output.decode('utf-8')))
        report_job_status = 'failed'
else:
    qstat_output = qstat_output.decode('utf-8')
    # logger.info('qstat call successful')
    report_job_status = parse_qstat_output(qstat_output, job_id)
finally:
    sys.stdout.write('{}\n'.format(report_job_status))
    # logger.info('<<< Done job id: {} / {}\n'.format(job_id, report_job_status))
    logging.shutdown()
    sys.exit(0)
