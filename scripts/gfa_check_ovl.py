#!/usr/bin/env python3

import re
import argparse as argp
import collections as col
import io
import statistics as stats


def parse_arguments():

    parser = argp.ArgumentParser(prog="gfa_check_ovl.py")

    parser.add_argument(
        '--graph',
        '-g',
        type=str,
        dest='graph',
        help='Full path to GFA input file.'
    )
    parser.add_argument(
        '--gfa-spec',
        '-s',
        type=int,
        default=1,
        choices=[1],
        dest='spec',
        help='Version of GFA specification. Default: 1'
    )
    parser.add_argument(
        '--out-stats',
        '-os',
        type=str,
        default='/dev/stdout',
        dest='out_stats',
        help='Full path to output file for statistics output.'
    )
    parser.add_argument(
        '--out-discard',
        '-od',
        type=str,
        dest='out_discard',
        default='auto'
    )
    args = parser.parse_args()
    return args


def determine_wrong_links(segment_lengths, links):

    error_links = []
    error_lengths = []

    for segment_id, segment_links in links.items():
        segment_length = segment_lengths[segment_id]
        for (ln, link_length, other_segment) in segment_links:
            length_diff = segment_length - link_length
            if length_diff < 1:
                error_lengths.append(link_length)
                error_links.append((ln, segment_id, other_segment, link_length, length_diff))
    
    return error_links, error_lengths


def parse_gfa_v1(input_graph):

    orientation_map = {
        ('+', '+'): 'frw-frw',
        ('-', '-'): 'rev-rev',
        ('+', '-'): 'frw-rev',
        ('-', '+'): 'rev-frw',   
    }

    nucleotides = re.compile('[ACGT]+')

    overlap_lengths = []
    overlap_stats = col.Counter()
    segment_lengths = dict()

    links = col.defaultdict(list)

    with open(input_graph, 'r') as gfa:
        for ln, line in enumerate(gfa, start=1):
            if not line[0] in ['S', 'L']:
                continue
            if line[0] == 'S':
                segment_len = None
                columns = line.strip().split('\t')
                try:
                    _, segment_id, seg_sequence = columns[:3]
                except ValueError:
                    print(ln)
                    print(line.strip())
                    raise
                remainder = columns[3:]
                if seg_sequence == '*':
                    for r in remainder:
                        if r.startswith('LN:'):
                            segment_len = int(r.split(':')[-1])
                else:
                    assert nucleotides.match(seg_sequence[:10].upper()) is not None, f'No nucleotide sequence: {seg_sequence}'
                    segment_len = len(seg_sequence)
                if segment_len is None:
                    raise ValueError(f'Could not determine segment length: {line.strip()}')
                segment_lengths[segment_id] = segment_len
            else:
                _, s1_id, s1_dir, s2_id, s2_dir, link_length = line.split('\t')[:6]
                link_length = int(link_length.strip().strip('M'))
                link_orientation = orientation_map[(s1_dir, s2_dir)]
                overlap_lengths.append(link_length)
                overlap_stats[link_orientation] += 1

                links[s1_id].append((ln, link_length, s2_id))
                links[s2_id].append((ln, link_length, s1_id))

    return segment_lengths, links, overlap_lengths, overlap_stats


def compute_length_statistics(length_values):

    length_values = sorted(length_values, reverse=True)
    total_length = sum(length_values)
    max_length = length_values[0]
    min_length = length_values[-1]
    median_length = int(round(stats.median(length_values), 0))
    mean_length = int(round(stats.mean(length_values), 0))
    n50_length = None
    cumsum = 0
    for l in length_values:
        cumsum += l
        if cumsum > total_length / 2:
            n50_length = l
            break

    return total_length, max_length, min_length, median_length, mean_length, n50_length


def compile_statistics(segment_lengths, overlap_lengths, overlap_stats, error_lengths):

    out_buffer = io.StringIO()

    total_length, max_length, min_length, median_length, mean_length, n50_length = compute_length_statistics(
        segment_lengths.values()
    )
    
    out_buffer.write(f'num_segments\t{len(segment_lengths)}\n')
    out_buffer.write(f'segments_total_length\t{total_length}\n')
    out_buffer.write(f'segments_min_length\t{min_length}\n')
    out_buffer.write(f'segments_max_length\t{max_length}\n')
    out_buffer.write(f'segments_median_length\t{median_length}\n')
    out_buffer.write(f'segments_mean_length\t{mean_length}\n')
    out_buffer.write(f'segments_N50_length\t{n50_length}\n')

    total_length, max_length, min_length, median_length, mean_length, n50_length = compute_length_statistics(
        overlap_lengths
    )

    out_buffer.write(f'num_links\t{len(overlap_lengths)}\n')
    out_buffer.write(f'links_total_length\t{total_length}\n')
    out_buffer.write(f'links_min_length\t{min_length}\n')
    out_buffer.write(f'links_max_length\t{max_length}\n')
    out_buffer.write(f'links_median_length\t{median_length}\n')
    out_buffer.write(f'links_mean_length\t{mean_length}\n')
    out_buffer.write(f'links_N50_length\t{n50_length}\n')

    out_buffer.write(f'links_num_errors\t{len(error_lengths)}\n')
    out_buffer.write(f'links_total_length_errors\t{sum(error_lengths)}\n')
    out_buffer.write(f'links_pct_len_errors\t{round(sum(error_lengths) / total_length * 100,3)}\n')
    for k in sorted(overlap_stats.keys()):
        out_buffer.write(f'links_orientation_{k}\t{overlap_stats[k]}\n')

    return out_buffer


def main():

    args = parse_arguments()
    if args.spec == 1:
        graph_parser = parse_gfa_v1
    else:
        raise NotImplementedError()

    segment_lengths, links, overlap_lengths, overlap_stats = graph_parser(args.graph)
    error_links, error_lengths = determine_wrong_links(segment_lengths, links)

    if args.out_discard == 'auto':
        out_discard = args.graph.replace('.gfa', '.discard.links')
    else:
        out_discard = args.out_discard

    with open(out_discard, 'w') as dump:
        _ = dump.write('# first line in GFA has index / line number 1\n')
        _ = dump.write('# GFA_LN - SEG1 - SEG2 - Link_Length - DiffLen_SEG1_LINK\n')
        _ = dump.write('\n'.join([f'{ln}\t{s1}\t{s2}\t{ll}\t{ld}' for (ln, s1, s2, ll, ld) in error_links]) + '\n')

    gfa_stats = compile_statistics(
        segment_lengths,
        overlap_lengths,
        overlap_stats,
        error_lengths
    )
    with open(args.out_stats, 'w') as dump:
        _ = dump.write(gfa_stats.getvalue())
    
    return


if __name__ == '__main__':
    main()
