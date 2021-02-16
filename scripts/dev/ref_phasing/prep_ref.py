#!/usr/bin/env python3

import os
import argparse
import re
import io


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta-in', '-i', dest='input', type=str)
    parser.add_argument('--fasta-out', '-o', dest='output', type=str)
    parser.add_argument('--chromosomes', '-c', dest='chrom', default='"^chr[0-9]+$"')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    chrom_match = re.compile(args.chrom.strip('"'))

    out_buffer = io.StringIO()

    collected_chroms = []

    collect = False
    with open(args.input, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                chrom = chrom_match.match(line.strip().strip('>'))
                if chrom is None:
                    collect = False
                    continue
                collected_chroms.append(chrom.group(0))
                out_buffer.write(line)
                collect = True
                continue
            elif collect:
                out_buffer.write(line)
            else:
                continue

    out_path = os.path.dirname(args.output)
    os.makedirs(out_path, exist_ok=True)
    with open(args.output, 'w') as dump:
        _ = dump.write(out_buffer.getvalue())

    print(sorted(collected_chroms))
    return 0

if __name__ == '__main__':
    main()