#!/usr/bin/env python

import os
import sys
import argparse
import re


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--outfile', '-o', type=str, dest='outfile')
    parser.add_argument('--at-least', '-a', type=str, dest='atleast')
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

    req_version = [int(v) for v in args.atleast.split('.')]

    version_pattern = re.compile('[0-9]+\\.[0-9]+(\\.[0-9]+)?')

    match_found = False

    with open(logfile, 'w') as log:
        _ = log.write('Minimum version required: {}\n'.format(args.atleast))
        for line in sys.stdin.readlines():
            _ = log.write('Processing line: {}\n'.format(line.strip()))
            mobj = version_pattern.search(line.strip())
            if mobj is not None:
                version_info = mobj.group(0)
                _ = log.write('Potential version info found: {}\n'.format(version_info))
                tool_version = [int(v) for v in version_info.split('.')]
                if all([x >= y for x, y in zip(req_version, tool_version)]):
                    _ = log.write('Minimum version matched...\n')
                    match_found = True
                    break
                else:
                    _ = log.write('Version info did not match...\n')

        if match_found:
            exit_code = 0
            with open(outfile, 'w') as touch:
                _ = touch.write('Version confirmed: {}\n'.format('.'.join([str(v) for v in tool_version])))
        else:
            exit_code = 1
            _ = log.write('No match found')

    return exit_code


if __name__ == '__main__':
    sys.exit(main())
