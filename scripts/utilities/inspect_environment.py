#!/usr/bin/env python

import os
import sys
import argparse
import subprocess as sp


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--outfile', '-o', type=str, dest='outfile')
    parser.add_argument('--logfile', '-l', type=str, dest='logfile')
    parser.add_argument('--export-conda-env', '-e', action='store_true', default=False, dest='export')

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

    conda_env = None

    with open(logfile, 'w') as log:

        _ = log.write('\n===== Accessible environment:\n')

        for k in env_vars:
            _ = log.write('{} - {}\n'.format(k, my_env[k]))
            if k == 'CONDA_PREFIX':
                conda_env = my_env[k]

        _ = log.write('\nDone\n')

        if args.export and conda_env is None:
            _ = logfile.write('\nERROR: cannot export CONDA env, no prefix path found in environment (see above)\n')
        elif args.export:
            _ = log.write('\n===== Export of active CONDA environment\n\n')

            try:
                out = sp.check_output('conda env export --prefix {}'.format(conda_env),
                                      stderr=sp.STDOUT,
                                      shell=True,
                                      env=None)
                out = out.decode('utf-8')
                _ = log.write(out + '\n\n')
            except sp.CalledProcessError as spe:
                _ = log.write('Exporting Conda env failed with code {}: {}\n'.format(spe.returncode, spe.output))
        else:
            pass

    with open(outfile, 'w') as touch:
        _ = touch.write('ENV OK\n')

    return


if __name__ == '__main__':
    main()
    sys.exit(0)
