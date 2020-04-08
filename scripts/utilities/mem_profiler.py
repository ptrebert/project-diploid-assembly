#!/usr/bin/env python3

import os
import time
import psutil


workdir = os.getcwd()
logfile = os.path.join(workdir, 'memlog.txt')

bytes_to_gigabytes = 1024 ** 3

time_limit = 86400

with open(logfile, 'w') as foo:
    pass

sleep_time = 0

with open(logfile, 'a') as log:
    header = '\t'.join(['#time', 'threads', 'load', 'mem_tot', 'mem_free', 'swap_tot', 'swap_free'])
    _ = log.write(header + '\n')
    while sleep_time < time_limit:
        now = str(time.ctime()).replace(' ', '_')
        threads = str(psutil.cpu_count(logical=True))
        pct_cpu = str(round(psutil.cpu_percent(), 2))
        mem = psutil.virtual_memory()
        mem_tot = str(round(mem.total / bytes_to_gigabytes, 2))
        mem_av = str(round(mem.available / bytes_to_gigabytes, 2))
        swap = psutil.swap_memory()
        swap_tot = str(round(swap.total / bytes_to_gigabytes, 2))
        swap_free = str(round(swap.free / bytes_to_gigabytes, 2))

        logline = '\t'.join([now, threads, pct_cpu, mem_tot, mem_av, swap_tot, swap_free])
        _ = log.write(logline + '\n')

        sleep_time += 60
        time.sleep(60)
