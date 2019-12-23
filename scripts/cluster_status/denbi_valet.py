#!/usr/bin/env python3

import os
import sys
import re
import logging
from logging.handlers import RotatingFileHandler
import subprocess as sp


log_msg_format = "%(asctime)s | %(levelname)s | %(message)s"
cluster_status_log = '/beeond/projects/diploid-assembly/log/cluster_status.log'
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


def parse_qstat_output(job_info):
    """
    :param job_info:
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
        'C': 'done',
        'E': 'done',
        'H': 'running',
        'M': 'running',
        'Q': 'running',
        'R': 'running',
        'S': 'running',
        'T': 'running',
        'W': 'running'
    }

    log_info, job_status, exit_code = '\n', None, -1
    for line in job_info.split('\n'):
        line = line.strip()
        if not line:
            continue
        if 'Job Id:' in line:
            logger.info('Processing status for job: {}'.format(line.split()[-1]))
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
                if exit_code != 0 and job_status == 'done':
                    job_status = 'failed'
                elif exit_code == 0 and job_status == 'done':
                    job_status = 'success'
                else:
                    if exit_code == 271:
                        logger.error('Job has exit code 271 - killed by PBS scheduler for exceeding resources '
                                     'or walltime limit.')
                    else:
                        logger.error('Job has exit code {}, but not status done: {}'.format(exit_code, job_status))
                    raise RuntimeError('Job has exit code but not status done')
        else:
            continue
    logger.info(log_info)
    if job_status == 'running':
        if exit_code != -1:
            logger.error('Job determined as running, but has exit code: {}'.format(exit_code))
            raise RuntimeError('Running job has exit code')
    return job_status


job_id = sys.argv[1]

logger.info('>>> Received input job id: {}'.format(job_id))

mobj = re.match('[0-9]+', job_id)

if mobj is None:
    logger.error('Could not determine numeric part of job ID - exiting...')
    raise ValueError('Malformed job id')

num_job_id = mobj.group(0)
report_job_status = 'failed'
try:
    qstat_output = sp.check_output(
        'qstat -Gf {}'.format(num_job_id),
        stderr=sp.STDOUT,
        shell=True,
        timeout=60
    )
except sp.CalledProcessError as cpe:
    logger.error('qstat call failed: {} / {}'.format(cpe.returncode, cpe.output.decode('utf-8')))
    if cpe.returncode == 153:
        # job id no longer known to qstat, assume it was fine...
        report_job_status = 'success'
    else:
        report_job_status = 'failed'
else:
    qstat_output = qstat_output.decode('utf-8')
    logger.info('qstat call successful')
    report_job_status = parse_qstat_output(qstat_output)
finally:
    sys.stdout.write('{}\n'.format(report_job_status))
    logger.info('<<< Done job id: {} / {}\n'.format(job_id, report_job_status))
    logging.shutdown()
    sys.exit(0)
