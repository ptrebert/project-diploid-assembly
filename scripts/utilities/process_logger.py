#!/usr/bin/env python3

import time

import psutil

attributes = [
    'cmdline',
    'cpu_percent',
    'create_time',
    'cwd',
    'exe',
    'memory_info',
    'name',
    'pid',
    'ppid',
    'status',
    'threads',
    'terminal',
    'username',
    'uids'
]

system_exe = [
    '/bin',
    '/usr/bin',
    '/lib',
    '/usr/lib',
    '/usr/sbin'
]

special_processes = [
    'sd-pam',
    'ssh-agent'
]

whitelist = [
    '/smk_env/'
]


def main():

    LOGFILE = '/home/ebertp/process.log'

    with open(LOGFILE, 'w') as logfile:
        pass

    while 1:
        cache = dict()
        suspects = []
        for process in psutil.process_iter(attrs=attributes, ad_value='N/A'):
            cache[process.info['pid']] = process.info['exe'], process.info['cmdline']
            if process.info['username'] != 'pebert':
                continue
            if any([sp in process.info['cmdline'][0] for sp in special_processes]):
                continue
            if any([process.info['exe'].startswith(se) for se in system_exe]):
                continue
            if any([wl in process.info['exe'] for wl in whitelist]):
                continue
            suspects.append(process.info)

        with open(LOGFILE, 'a') as logfile:
            for p_info in suspects:
                _ = logfile.write('PARENT: {} / {}\n'.format(*cache[p_info['ppid']]))
                _ = logfile.write('OFFENDER\n')
                block = '\n'.join(['{}\t{}'.format(k, p_info[k]) for k in sorted(p_info.keys())])
                _ = logfile.write(block + '\n')
                _ = logfile.write('========\n')

        time.sleep(30)
    return


if __name__ == '__main__':
    main()
