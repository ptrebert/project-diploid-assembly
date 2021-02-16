#!/usr/bin/env python3

import os
import io
import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('tags', type=str)
    args = parser.parse_args()
    return args


def check_haplotags(file_path):

    names = ['readname', 'haplotype', 'phaseset', 'chromosome']
    df = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=names)
    
    hap_counts = df['haplotype'].value_counts()
    total = hap_counts.sum()
   
    print('--- percent tagged')

    for hap, count in hap_counts.items():
        print(hap, round(count / total * 100, 2))


    h1_reads = set(df.loc[df['haplotype'] == 'H1', 'readname'].values)
    h2_reads = set(df.loc[df['haplotype'] == 'H2', 'readname'].values)
    untagged = set(df.loc[df['haplotype'] == 'none', 'readname'].values)

    print('--- intersect')

    print('h1 v h2 ', len(h1_reads.intersection(h2_reads)))
    print('h1 v h0 ', len(h1_reads.intersection(untagged)))
    print('h2 v h0 ', len(h2_reads.intersection(untagged)))

    return


def main():
    args = parse_args()
    check_haplotags(args.tags)

    return 0


if __name__ == '__main__':
    main()
