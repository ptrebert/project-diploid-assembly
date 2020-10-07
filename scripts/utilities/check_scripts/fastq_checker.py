#!/usr/bin/env python3

import os
import io
import argparse

import dnaio

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq1', type=str)
    parser.add_argument('fastq2', type=str)
    args = parser.parse_args()
    return args


def collect_read_names(fastq_path):

    reads = []

    with dnaio.open(fastq_path) as fastx:
        for record in fastx:
            reads.append(record.name)

    total_reads = len(reads)
    reads = set(reads)

    if not total_reads == len(reads):
        print('error: read duplicates {}: {} out of {}'.format(fastq_path, total_reads - len(reads), total_reads))

    return reads


def main():
    args = parse_args()
    fq1_reads = collect_read_names(args.fastq1)
    fq2_reads = collect_read_names(args.fastq2)

    intersect = fq1_reads.intersection(fq2_reads)
    if len(intersect) > 0:
        print('error: read sets not disjoined: {} out of {} / {}'.format(len(intersect), len(fq1_reads), len(fq2_reads)))

    print('fq1 reads {}'.format(len(fq1_reads)))
    print('fq2 reads {}'.format(len(fq2_reads)))


    return 0


if __name__ == '__main__':
    main()
