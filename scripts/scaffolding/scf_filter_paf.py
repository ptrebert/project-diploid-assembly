#!/usr/bin/env python3

import sys as sys
import os as os
import argparse as argp
import traceback as trb

import pandas as pd

PAF_HEADER = [
    'query_name',
    'query_length',
    'aln_qstart',
    'aln_qend',
    'aln_qstrand',
    'target_name',
    'target_length',
    'aln_tstart',
    'aln_tend',
    'aln_num_match',
    'aln_block_length',
    'aln_mapq',
    'aln_type',
    'aln_num_minimizer',
    'aln_chaining_score',
    'aln_chaining_score2',
    'query_region_repseeds'
]

PAF_OUTPUT_HEADER = PAF_HEADER[:12]

NUMERIC_HINTS = [
    'length',
    'start',
    'end',
    'mapq',
    'match'
]

PAF_COLUMN_DTYPES = dict((c, int) if any([x in c for x in NUMERIC_HINTS]) else (c, str) for c in PAF_HEADER)


def parse_command_line():

    parser = argp.ArgumentParser()
    parser.add_argument(
        '--paf',
        '-p',
        dest='paf',
        type=str,
        help='Full path to PAF file.'
    )
    parser.add_argument(
        '--output',
        '-o',
        dest='output',
        type=str,
        help='Full path to output file.'
    )
    parser.add_argument(
        '--min-matched-residues',
        '-m',
        type=int,
        default=500,
        dest='min_match',
        help='Matched residues threshold (lower limit)'
    )
    parser.add_argument(
        '--min-mapq-threshold',
        '-q',
        type=int,
        default=60,
        dest='min_mapq',
        help='MAPQ threshold (lower limit)'
    )
    parser.add_argument(
        '--tip-boundary',
        '-b',
        type=int,
        default=5,
        dest='tip_boundary',
        help='Tip boundary as percent of query/target length.'
    )
    parser.add_argument(
        '--containment',
        '-c',
        type=int,
        default=90,
        dest='containment',
        help='Containment threshold for query aligned length (upper limit).'
    )
    args = parser.parse_args()
    return args


def check_mapq(lower_threshold, alignments):
    return alignments.loc[alignments['aln_mapq'] >= lower_threshold, :]


def check_matched_residues(lower_threshold, alignments):
    return alignments.loc[alignments['aln_num_match'] >= lower_threshold, :]


def check_contained_alignment(upper_threshold, alignments):
    
    query_pct_size = (alignments['query_length'] * upper_threshold).round(0)
    query_aligned_length = alignments['aln_qend'] - alignments['aln_qstart']
    assert (query_aligned_length > 0).all(), 'Q aligned length wrong contains 0'
    selector = query_aligned_length < query_pct_size
    return alignments.loc[selector, :]


def check_connecting_alignment(boundary, alignments):

    alignments['reason'] = 'none'

    query_pct_size = (alignments['query_length'] * boundary).round(0)
    query_left = alignments['aln_qstart'] < query_pct_size
    query_right = alignments['aln_qend'] > alignments['query_length'] - query_pct_size
    query_not_connecting = ~(query_left | query_right)
    alignments.loc[query_not_connecting, 'reason'] = 'Qconnect'

    target_pct_size = (alignments['target_length'] * boundary).round(0)
    target_left = alignments['aln_tstart'] < target_pct_size
    target_right = alignments['aln_tend'] > alignments['target_length'] - target_pct_size
    target_not_connecting = ~(target_left | target_right)

    # select where query is connecting, but target is not
    alignments.loc[~query_not_connecting & target_not_connecting, 'reason'] = 'Tconnect'

    reasons = alignments['reason'].value_counts()

    alignments = alignments.loc[alignments['reason'] == 'none', ].copy()
    alignments.drop('reason', axis=1, inplace=True)

    return alignments, reasons


def check_singleton(alignments):

    query_counts = alignments['query_name'].value_counts()
    singletons = set(query_counts.index[query_counts == 1])

    return alignments.loc[~alignments['query_name'].isin(singletons), :]


def check_local_alignment(alignments):
    """
    Read has only (but several) alignments on
    a single contig.
    """
    discard = set()
    for query_name, subset in alignments.groupby(['query_name']):
        if subset['target_name'].nunique() == 1:
            discard.add(query_name)
    return alignments.loc[~alignments['query_name'].isin(discard), :]


def check_multi_alignment(alignments):
    """
    Read has several alignments to one contig,
    select best one in terms of res matches
    """
    select_indices = set()
    for (read_name, contig_name), sub_aln in alignments.groupby(['query_name', 'target_name']):
        if sub_aln.shape[0] > 1:
            selector = sub_aln['aln_num_match'] == sub_aln['aln_num_match'].max()
            select_index = sub_aln.index[selector]
            if select_index.size != 1:
                # then select shorter block length
                selector &= sub_aln['aln_block_length'] == sub_aln['aln_block_length'].min()
                select_index = sub_aln.index[selector]
                if select_index.size != 1:
                    err_msg = 'Several alignments with identical match score and block length: {} to {}\n'.format(read_name, contig_name)
                    err_msg += 'Max match score: {}\n'.format(sub_aln['aln_num_match'].max())
                    err_msg += 'Selected alignments: {}\n'.format(select_index.size)
                    err_msg += 'Total subset size: {}\n'.format(sub_aln.shape)
                    raise ValueError(err_msg)
            select_indices = select_indices.union(set(select_index))
        else:
            select_indices = select_indices.union(set(sub_aln.index))
    return alignments.loc[alignments.index.isin(select_indices), :]



def main():

    args = parse_command_line()

    df = pd.read_csv(
        args.paf,
        sep='\t',
        header=None,
        names=PAF_HEADER,
        usecols=PAF_OUTPUT_HEADER,
        dtype=PAF_COLUMN_DTYPES
    )

    stats = [('total_alignments', df.shape[0])]
    last_num_records = df.shape[0]

    df = check_mapq(args.min_mapq, df)
    discarded = last_num_records - df.shape[0]
    stats.append(('discard_mapq', discarded))
    last_num_records = df.shape[0]

    df = check_matched_residues(args.min_match, df)
    discarded = last_num_records - df.shape[0]
    stats.append(('discard_Nmatches', discarded))
    last_num_records = df.shape[0]

    df = check_contained_alignment(args.containment / 100, df)
    discarded = last_num_records - df.shape[0]
    stats.append(('discard_containment', discarded))
    last_num_records = df.shape[0]

    # here, we need a copy because columns are added to the df
    df, reasons = check_connecting_alignment(args.tip_boundary / 100, df.copy())
    stats.append(('discard_Qconnect', reasons['Qconnect']))
    stats.append(('discard_Tconnect', reasons['Tconnect']))
    last_num_records = df.shape[0]

    df = check_singleton(df)
    discarded = last_num_records - df.shape[0]
    stats.append(('discard_singleton', discarded))
    last_num_records = df.shape[0]

    df = check_local_alignment(df)
    discarded = last_num_records - df.shape[0]
    stats.append(('discard_local', discarded))
    last_num_records = df.shape[0]

    df = check_multi_alignment(df)
    discarded = last_num_records - df.shape[0]
    stats.append(('discard_multi', discarded))

    stats.append(('filtered_alignments', df.shape[0]))

    output_folder = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(output_folder, exist_ok=True)

    df.to_csv(
        args.output,
        sep='\t',
        header=True,
        index=False
    )

    stat_output = args.output.rsplit('.', 1)[0]
    if stat_output.endswith('.paf'):
        # compressed output
        stat_output = stat_output.rsplit('.', 1)[0]
    stat_output += '.stats'
    with open(stat_output, 'w') as dump:
        _ = dump.write('\n'.join(['{}\t{}'.format(k, v) for k, v in stats]) + '\n')

    return


if __name__ == '__main__':
    try:
        main()
    except Exception:
        trb.print_exc(file=sys.stderr)
        raise
    sys.exit(0)
