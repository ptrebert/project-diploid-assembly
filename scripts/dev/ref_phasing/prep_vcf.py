#!/usr/bin/env python3

import os
import argparse
import re
import io

import pysam

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf-in', '-i', dest='input', type=str)
    parser.add_argument('--out-pattern', '-o', dest='output', type=str)
    parser.add_argument('--chromosomes', '-c', dest='chrom', default='"^chr[0-9]+$"')
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    chrom_match = re.compile(args.chrom.strip('"'))

    with pysam.VariantFile(args.input) as vcf:
        collected_chroms = [c for c in vcf.header.contigs if chrom_match(c) is not None]

    call = 'bcftools view --regions {} --output-type v --output-file {} {}'
    for c in collected_chroms:
        out_file = args.output.strip('"').format(c)
        out_path = os.path.dirname(out_file)
        os.makedirs(out_path, exist_ok=True)
        tmp = call.format(c, out_file, args.input)
        try:
            out = subprocess.check_output(tmp, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')
        except subprocess.CalledProcessError as spe:
            raise RuntimeError(spe.output.decode('utf-8'))

    print(sorted(collected_chroms))
    return 0

if __name__ == '__main__':
    main()