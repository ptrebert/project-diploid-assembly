#!/usr/bin/env python

import os
import sys
import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--outfile', '-o', type=str, dest='outfile')
    parser.add_argument('--logfile', '-l', type=str, dest='logfile')

    args = parser.parse_args()

    outfile = args.outfile
    logfile = args.logfile

    try:
        os.makedirs(os.path.dirname(os.path.abspath(outfile)), exist_ok=True)
        os.makedirs(os.path.dirname(os.path.abspath(logfile)), exist_ok=True)
    except TypeError:
        # since Conda environments (or the Singularity module on Hilbert)
        # only support Python2 (...), exist_ok may cause an exception
        # Ignore that and hope that Snakemake creates everything...
        pass

    my_env = dict(os.environ)

    env_vars = sorted(my_env.keys())

    with open(logfile, 'w') as log:

        _ = log.write('\nAccessible environment:\n')

        for k in env_vars:
            _ = log.write('{} - {}\n'.format(k, my_env[k]))

        _ = log.write('\nDone\n')

    with open(outfile, 'w') as touch:
        _ = touch.write('ENV OK\n')

    return


if __name__ == '__main__':
    main()
    sys.exit(0)
