#!/usr/bin/env python3

import sys
import fileinput
import argparse as argp


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        '--lib-name',
        '-l',
        type=str,
        dest='libname',
        required=True
    )
    parser.add_argument(
        '--separator',
        '-s',
        type=str,
        dest='separator',
        default='_'
    )
    parser.add_argument(
        '--keep-ambig',
        '-k',
        action='store_true',
        default=False,
        dest='keep_ambig'
    )
    args = parser.parse_args()
    return args


def main():

    args = parse_command_line()
    libname = args.libname.strip()

    process_header = lambda line: line.strip() + '_' + libname + '\n'
    if args.keep_ambig:
        process_seq = lambda line: line
    else:
        process_seq = lambda line: line.strip().strip('N') + '\n'
    
    for line in fileinput.input([]):
        if not line.strip():
            continue
        if line.startswith('>'):
            sys.stdout.write(process_header(line))
        else:
            sys.stdout.write(process_seq(line))
    return


if __name__ == '__main__':
    main()
