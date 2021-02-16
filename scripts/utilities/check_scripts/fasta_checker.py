#!/usr/bin/env python3

import os
import io
import argparse
import hashlib

import dnaio

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta1', type=str)
    parser.add_argument('fasta2', type=str)
    args = parser.parse_args()
    return args


def collect_contig_hashes(fasta_path):

    contigs = dict()
    duplicates = list()

    with dnaio.open(fasta_path) as fastx:
        for record in fastx:
            contig_name = record.name
            seq = record.sequence
            seq_hash = hashlib.blake2b(seq.encode('ascii')).digest()
            if seq_hash in contigs:
                # happens for flye assembler
                duplicates.append(contig_name)
            contigs[seq_hash] = contig_name, len(seq)

    print('assembled contigs {}'.format(len(contigs)))
    if duplicates:
        print('duplicate contig sequences: {}'.format(sorted(duplicates)))

    return contigs


def main():
    args = parse_args()
    fasta1 = collect_contig_hashes(args.fasta1)
    fasta2 = collect_contig_hashes(args.fasta2)

    total1 = sum([x[1] for x in fasta1.values()])
    total2 = sum([x[1] for x in fasta2.values()])

    intersect = set(fasta1.keys()).intersection(fasta2.keys())
    if len(intersect) > 0:
        dups1 = [fasta1[h] for h in intersect]
        dups2 = [fasta2[h] for h in intersect]

        cluster1 = sum([x[1] for x in dups1])
        cluster2 = sum([x[1] for x in dups2])

        frac_h1 = round(len(dups1) / len(fasta1) * 100, 2)
        frac_h2 = round(len(dups2) / len(fasta2) * 100, 2)

        pct_bp_h1 = round(cluster1 / total1 * 100, 2)
        pct_bp_h2 = round(cluster2 / total2 * 100, 2)

        print('HOM frac h1 #ctg', frac_h1)
        print('HOM frac h1 pct. bp', pct_bp_h1)

        print('HOM frac h2 #ctg', frac_h2)
        print('HOM frac h2 pct. bp', pct_bp_h2)

    return 0


if __name__ == '__main__':
    main()
