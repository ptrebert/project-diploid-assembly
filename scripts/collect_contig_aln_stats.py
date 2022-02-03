#!/usr/bin/env python

import os
import sys
import pathlib as pl
import argparse
import csv
import operator as op
import collections as col
import re

import pandas as pd

def parse_command_line():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--contig-alignments',
        '-ca',
        type=lambda x: pl.Path(x).resolve(),
        dest='alignments',
        help='Contig to reference alignment as BED file'
    )
    parser.add_argument(
        '--reference-chromosomes',
        '-rc',
        type=lambda x: pl.Path(x).resolve(),
        dest='refchroms',
        help='Reference chromosome names and lengths as 2-column TSV table or FASTA index (fai) file.'
    )
    parser.add_argument(
        '--contig-names',
        '-cn',
        type=lambda x: pl.Path(x).resolve(),
        dest='contigs',
        help='Contig names and lengths as 2-column TSV table or FASTA index (fai) file.'
    )
    parser.add_argument(
        '--min-mapq',
        '-m',
        type=int,
        default=0,
        dest='minmapq'
    )
    parser.add_argument(
        '--min-chrom-coverage',
        '-mcc',
        type=int,
        default=0,
        choices=range(0, 101),
        metavar='[0...100]',
        dest='min_chrom_coverage',
        help='Specify the minimum (alignment) coverage per chromosome as percent of chromosome length '
                'covered that has to be attained for a single cluster. Otherwise, this script will '
                'abort as it is assumed that there is at least one large assembly error for this '
                'cluster/chromosome. Default: 0 pct. (0...100)'
    )
    parser.add_argument(
        '--ignore-min-chrom-cov',
        '-ignore',
        action='store_true',
        default=False,
        dest='ignore_coverage_error',
        help='If set, ignore reference chromosome alignment coverages below the minimal threshold, i.e. '
             'do not fail and do not raise an error. If not set [the default], dump the alignment table '
             'with suffix ".error" and raise an exception.'
    )
    parser.add_argument(
        '--chrom-coverage-select',
        '-ccs',
        type=str,
        default='^chr[0-9X]+$',
        dest='chrom_coverage_select',
        help='Specify regular expression to limit checking the minimum coverage threshold only for '
                'matching chromosomes. Default: chr[0-9X]+ (i.e., autosomes and chrX)'
    )
    parser.add_argument(
        '--combine-sequences',
        '-comb',
        nargs='*',
        default=None,
        dest='combine',
        help='Combine these sequences from the reference into a single statistic.'
    )
    parser.add_argument(
        '--contig-groups',
        '-cg',
        action='store_true',
        default=False,
        dest='groups',
        help='Assume individual contigs are organized in groups (scaffolds, cluster). Default: False'
    )
    parser.add_argument(
        '--group-id-position',
        '-gip',
        type=int,
        default=None,
        nargs='*',
        dest='groupid_pos',
        help='If contigs are grouped, the group id is composed of these parts '
             'in the BED name field (0-based indexing).'
    )
    parser.add_argument(
        '--group-id-split',
        '-gis',
        type=str,
        dest='split_char',
        default='_',
        help='Use this character for splitting the BED name column to extract the group id. Default: "_"'
    )
    parser.add_argument(
        '--output',
        '-o',
        type=lambda x: pl.Path(x),
        dest='output',
        help='Path to TSV output table.'
    )
    args = parser.parse_args()
    return args


def build_name_extract_function(args):

    def get_name(x):
        return x
    if args.groups:
        groupid_pos = args.groupid_pos
        if groupid_pos is None:
            raise ValueError('Contigs are grouped but no positional index for group ID specified')
        elif isinstance(groupid_pos, int):
            groupid_pos = [groupid_pos]
        else:
            pass
        get_parts = op.itemgetter(*tuple(groupid_pos))

        def get_name(x):
            return ''.join(get_parts(x.split(args.split_char)))
    return get_name


def read_chromosome_sizes(fpath, combine_chroms, name_extractor):

    chrom_sizes = col.Counter()
    genome_size = 0
    chrom_names = {}
    with open(fpath, 'r') as table:
        for line in table:
            columns = line.split()
            seq_name = columns[0]
            seq_extracted_name = name_extractor(columns[0])
            seq_length = int(columns[1])

            chrom_sizes[seq_name] += seq_length
            genome_size += seq_length
            chrom_names[columns[0]] = seq_extracted_name
            if seq_name != seq_extracted_name:
                chrom_sizes[seq_extracted_name] += seq_length

    if combine_chroms is not None:
        combined_length = sum(chrom_sizes[c] for c in combine_chroms)
        if combined_length > 0:
            chrom_sizes[''.join(sorted(combine_chroms))] += combined_length
            for c in combine_chroms:
                chrom_names[c] = ''.join(sorted(combine_chroms))
    chrom_sizes['genome'] = genome_size
    return chrom_sizes, chrom_names


def read_contig_alignment_table(fpath, ref_seqs, assm_seqs, minmapq):
    """
    Interval merging based on nifty trick found in
    https://stackoverflow.com/a/65282946
    """

    strand_mapping = {
        '+': 1,
        '-': -1,
        '.': 0
    }

    columns = [
        'ref_seq',
        'start',
        'end',
        'assm_seq',
        'mapq',
        'strand'
    ]

    df = pd.read_csv(
        fpath,
        sep='\t',
        header=None,
        names=columns
    )
    df = df.loc[df['mapq'] >= minmapq, :].copy()

    # if user specified, this would replace, e.g., chrX and chrY with chrXchrY
    df['ref_seq'].replace(ref_seqs, inplace=True)

    # if contigs are grouped (i.e. by cluster ID),
    # this replaces the full name by the group [cluster] ID
    df['assm_seq'].replace(assm_seqs, inplace=True)

    df['strand'].replace(strand_mapping, inplace=True)

    df.sort_values(['ref_seq', 'start', 'end'], inplace=True, ascending=True)

    aln_lengths = col.Counter()
    for (ref_seq, assm_seq), alignments in df.groupby(['ref_seq', 'assm_seq']):
        intervals = pd.arrays.IntervalArray.from_tuples(
            [(row.start, row.end) for row in alignments.itertuples(index=False)],
        )
        if intervals.is_non_overlapping_monotonic:
            aln_length = sum([iv.length for iv in intervals])
            aln_lengths[(ref_seq, assm_seq)] += aln_length
            continue
        merge_intervals = alignments[['start', 'end']].copy()
        merge_intervals['interval_id'] = (merge_intervals['start'] > merge_intervals['end'].shift().cummax()).cumsum()
        merge_intervals = merge_intervals.groupby("interval_id").agg({"start":"min", "end": "max"})
        
        aln_length = (merge_intervals['end'] - merge_intervals['start']).sum()
        aln_lengths[(ref_seq, assm_seq)] += aln_length

    return aln_lengths


def collect_row_alignment_stats(ref_seq, aln_lengths, ref_seq_sizes, assm_seq_sizes, min_chrom_cov, chrom_cov_select):

    if ref_seq == 'genome':
        ref_counts = col.Counter()
        for (ref_seq_name, assm_seq_name), aln_bp in aln_lengths.items():
            ref_counts[assm_seq_name] += aln_bp
        ref_counts['genome'] = sum(ref_counts.values())
        top_keys = ['top1_alignment']
    else:
        ref_counts = col.Counter({k[1]: v for k, v in aln_lengths.items() if k[0] == ref_seq})
        top_keys = ['top1_alignment', 'top2_alignment', 'top3_alignment']

    row_records = dict()
    ref_counts = ref_counts.most_common()

    error_record = ''

    chrom_matcher = re.compile(chrom_cov_select.strip('"'))

    for k, (assm_seq, aligned_bases) in zip(top_keys, ref_counts[:3]):
        assm_pct = str(round(aligned_bases / assm_seq_sizes[assm_seq] * 100, 2))
        ref_pct = round(aligned_bases / ref_seq_sizes[ref_seq] * 100, 2)
        if k == 'top1_alignment':
            # change here: to make (i) debugging easier, or to (ii) ignore the chrom coverage error,
            # dump the table either (i) to an alternate location (change file extension) and raise
            # the below error, or (ii) dump to specified output location and print as warning
            if chrom_matcher.match(ref_seq) is not None and ref_pct < min_chrom_cov:
                error_record += f'Top 1 alignment from assembly {assm_seq} below coverage threshold '
                error_record += f'for reference chromosome {ref_seq}: {ref_pct}% (< {min_chrom_cov}%)'

        item = '|'.join([assm_seq, 'bp:' + str(aligned_bases), 'asm:' + assm_pct, 'ref:' + str(ref_pct)])
        row_records[k] = item

    ref_counts = sorted(ref_counts, key=lambda x: x[0])
    for assm_seq, aligned_bases in ref_counts:
        assm_pct = str(round(aligned_bases / assm_seq_sizes[assm_seq] * 100, 2))
        ref_pct = str(round(aligned_bases / ref_seq_sizes[ref_seq] * 100, 2))
        item = '|'.join(['bp:' + str(aligned_bases), 'asm:' + assm_pct, 'ref:' + ref_pct])
        row_records[assm_seq] = item
    return row_records, error_record


def create_output_table(aln_lengths, ref_seq_sizes, assm_seq_sizes, min_chrom_cov, chrom_cov_select):

    aligned_contigs = sorted(set(k[1] for k in aln_lengths.keys()))
    aligned_refs = sorted(set(k[0] for k in aln_lengths.keys()), key=lambda x: ref_seq_sizes[x], reverse=True)

    table_rows = []
    error_records = []
    for ref_seq in ref_seq_sizes.keys():
        row_record, error_record = collect_row_alignment_stats(
            ref_seq,
            aln_lengths,
            ref_seq_sizes,
            assm_seq_sizes,
            min_chrom_cov,
            chrom_cov_select
        )
        row_record['ref_seq'] = ref_seq
        row_record['ref_length'] = ref_seq_sizes[ref_seq]
        table_rows.append(row_record)
        if error_record:
            error_records.append(error_record)

    out_header = [
        'ref_seq',
        'ref_length',
        'top1_alignment',
        'top2_alignment',
        'top3_alignment',
    ]
    out_header.extend(aligned_contigs)

    return table_rows, out_header, error_records


def main():
    args = parse_command_line()
    name_extractor = build_name_extract_function(args)
    ref_chroms_sizes, ref_chrom_names = read_chromosome_sizes(
        args.refchroms,
        args.combine,
        lambda x: x
    )
    assm_chrom_sizes, assm_chrom_names = read_chromosome_sizes(
        args.contigs,
        [],
        name_extractor
    )
    aln_lengths = read_contig_alignment_table(
        args.alignments,
        ref_chrom_names,
        assm_chrom_names,
        args.minmapq
    )
    table_rows, header, error_records = create_output_table(
        aln_lengths,
        ref_chroms_sizes,
        assm_chrom_sizes,
        args.min_chrom_coverage,
        args.chrom_coverage_select
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    err_msg = '{} - alignment coverage [{}] below minimal threshold for the following reference chromosomes:\n'
    err_msg += '\n'.join(error_records) + '\n'
    if args.ignore_coverage_error or not error_records:
        table_file = args.output
        err_msg = err_msg.format('WARNING', args.alignments.name)
    else:
        table_file = args.output.with_suffix('.tsv.error')
        err_msg = err_msg.format('ERROR', args.alignments.name)

    with open(table_file, 'w', newline='') as table:
        _ = table.write('## ALIGN RECORDS: bp - base pair // asm: pct. aligned assembly sequence // ref: pct. aligned reference sequence\n')
        _ = table.write('#')
        writer = csv.DictWriter(table, fieldnames=header, quoting=csv.QUOTE_NONE,
                                extrasaction='ignore', restval='NA', delimiter='\t')
        writer.writeheader()
        writer.writerows(table_rows)

    if error_records and table_file == args.output:
        sys.stderr.write(err_msg)
    elif error_records:
        raise ValueError(err_msg)
    else:
        pass

    return 0


if __name__ == '__main__':
    main()
