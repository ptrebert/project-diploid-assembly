#!/usr/bin/env python

import os
import argparse
import csv
import operator as op
import collections as col


def parse_command_line():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--contig-alignments',
        '-ca',
        type=str,
        dest='alignments',
        help='Contig to reference alignment as BED file'
    )
    parser.add_argument(
        '--reference-chromosomes',
        '-rc',
        type=str,
        dest='refchroms',
        help='Reference chromosome names and lengths as 2-column TSV table or FASTA index (fai) file.'
    )
    parser.add_argument(
        '--contig-names',
        '-cn',
        type=str,
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
        type=str,
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

    strand_mapping = {
        '+': 1,
        '-': -1
    }
    aln_infos = col.defaultdict(list)
    aln_lengths = col.Counter()
    with open(fpath, 'r') as table:
        for line in table:
            columns = line.strip().split()
            ref_seq = columns[0]
            start = int(columns[1])
            end = int(columns[2])
            assm_seq = columns[3]  # BED name field
            mapq = int(columns[4])  # BED score field
            if mapq < minmapq:
                continue
            strand = strand_mapping.get(columns[5], 0)
            ref_seq_store = ref_seqs[ref_seq]
            assm_seq_store = assm_seqs[assm_seq]
            aligned_length = end - start
            aln_infos[(ref_seq_store, assm_seq_store)].append((aligned_length, mapq, strand))
            aln_lengths[(ref_seq_store, assm_seq_store)] += aligned_length
    return aln_infos, aln_lengths


def collect_row_alignment_stats(ref_seq, aln_lengths, ref_seq_sizes, assm_seq_sizes):

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

    for k, (assm_seq, aligned_bases) in zip(top_keys, ref_counts[:3]):
        assm_pct = str(round(aligned_bases / assm_seq_sizes[assm_seq] * 100, 2))
        ref_pct = str(round(aligned_bases / ref_seq_sizes[ref_seq] * 100, 2))
        item = '|'.join([assm_seq, str(aligned_bases), assm_pct, ref_pct])
        row_records[k] = item

    ref_counts = sorted(ref_counts, key=lambda x: x[0])
    for assm_seq, aligned_bases in ref_counts:
        assm_pct = str(round(aligned_bases / assm_seq_sizes[assm_seq] * 100, 2))
        ref_pct = str(round(aligned_bases / ref_seq_sizes[ref_seq] * 100, 2))
        item = '|'.join([str(aligned_bases), assm_pct, ref_pct])
        row_records[assm_seq] = item
    return row_records


def create_output_table(aln_info, aln_lengths, ref_seq_sizes, assm_seq_sizes):

    aligned_contigs = sorted(set(k[1] for k in aln_lengths.keys()))
    aligned_refs = sorted(set(k[0] for k in aln_lengths.keys()), key=lambda x: ref_seq_sizes[x], reverse=True)

    table_rows = []
    for ref_seq in aligned_refs:
        row_record = collect_row_alignment_stats(
            ref_seq,
            aln_lengths,
            ref_seq_sizes,
            assm_seq_sizes
        )
        row_record['ref_seq'] = ref_seq
        row_record['ref_length'] = ref_seq_sizes[ref_seq]
        table_rows.append(row_record)

    row_record = collect_row_alignment_stats(
        'genome',
        aln_lengths,
        ref_seq_sizes,
        assm_seq_sizes
    )
    row_record['ref_seq'] = 'genome'
    row_record['ref_length'] = ref_seq_sizes['genome']
    table_rows.append(row_record)

    out_header = [
        'ref_seq',
        'ref_length',
        'top1_alignment',
        'top2_alignment',
        'top3_alignment',
    ]
    out_header.extend(aligned_contigs)

    return table_rows, out_header


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
    aln_info, aln_lengths = read_contig_alignment_table(
        args.alignments,
        ref_chrom_names,
        assm_chrom_names,
        args.minmapq
    )
    table_rows, header = create_output_table(
        aln_info,
        aln_lengths,
        ref_chroms_sizes,
        assm_chrom_sizes
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, 'w', newline='') as table:
        _ = table.write('#')
        writer = csv.DictWriter(table, fieldnames=header, quoting=csv.QUOTE_NONE,
                                extrasaction='ignore', restval='NA', delimiter='\t')
        writer.writeheader()
        writer.writerows(table_rows)
    return 0


if __name__ == '__main__':
    main()
