#!/usr/bin/env python3

import os as os
import sys as sys
import resource
import logging as log
import argparse as argp
import traceback as trb
import collections as col
import io as io

import xopen as xopen
import pandas as pd
import networkx as nx


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog="gfa_subset.py")
    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        default=False,
        dest="debug",
        help="Print status and error messages to STDOUT. Otherwise, only "
            "errors/warnings will be reported to STDERR.",
    )
    parser.add_argument(
        "--input-gfa",
        "-ig",
        type=str,
        required=True,
        dest="input_gfa",
        help="Full path GFA (can be compressed).",
    )
    parser.add_argument(
        "--input-table",
        "-it",
        type=str,
        required=True,
        dest="input_table",
        help="Full path to tig subset table."
    )
    parser.add_argument(
        "--output-gfa",
        "-og",
        type=str,
        required=True,
        dest="output_gfa",
        help="Full path to output (subset) GFA file."
    )
    parser.add_argument(
        "--output-table",
        "-ot",
        type=str,
        required=True,
        dest="output_table",
        help="Full path to CSV output table (to be loaded in Bandage)."
    )
    return parser.parse_args()


def buffer_graph(gfa_path):
    """
    """
    links = []
    node_map = []
    with xopen.xopen(gfa_path, 'r') as graph:
        for ln, line in enumerate(graph):
            line_type = line[0]
            if line_type == '#':
                continue
            elif line_type == 'S':
                node = line.split()[1]
                node_map.append((ln, node))
                continue
            elif line_type == 'A':
                continue
            elif line_type == 'L':
                _, node1, _, node2 = line.split()[:4]
                links.append((node1, node2))
                node_map.append((ln, node1))
                node_map.append((ln, node2))
            else:
                raise ValueError('Unexpected line type in GFA: {} / {}'.format(line_type, line[:30]))
    return links, node_map


def load_table(tig_table):
    """
    """
    header = [
        'aln_chrom',
        'aln_start',
        'aln_end',
        'tig_name',
        'mapq',
        'orientation',
        'reg_chrom',
        'reg_start',
        'reg_end',
        'reg_label',
        'reg_color',
        'overlap_bp'
    ]
    df = pd.read_csv(
        tig_table,
        sep='\t',
        header=None,
        names=header
    )

    select_tigs = set()
    annotations = col.defaultdict(list)

    tig_label_colors =  df.groupby(['tig_name', 'reg_label', 'reg_color'])['overlap_bp'].sum()
    for index_tig, overlaps in tig_label_colors.groupby('tig_name'):
        overlap_bp = overlaps.iloc[0]
        tig, label, color = overlaps.index[0]
        select_tigs.add(tig)
        annotations[tig].append((overlap_bp, label, color))

    return select_tigs, annotations


def determine_connected_tigs(edges, select_tigs):
    """
    """
    g = nx.Graph()
    g.add_edges_from(edges)

    num_cc = 0
    connected_tigs = set().union(select_tigs)
    for cc in nx.connected_components(g):
        num_cc += 1
        if len(select_tigs.intersection(cc)) == 0:
            continue
        connected_tigs = connected_tigs.union(cc)
    return connected_tigs, num_cc


def main(logger, cargs):
    """
    """
    logger.debug('Loading GFA file: {}'.format(cargs.input_gfa))
    links, node_map = buffer_graph(cargs.input_gfa)
    logger.debug('GFA loading complete')
    select_tigs, annotations = load_table(cargs.input_table)
    logger.debug('{} tigs selected as GFA subset'.format(len(select_tigs)))

    connected_tigs, num_cc = determine_connected_tigs(links, select_tigs)
    logger.debug('{} tigs selected as connected (total {} conn. comp.)'.format(len(connected_tigs), num_cc))

    keep_lines = set(t[0] for t in node_map if t[1] in connected_tigs)

    output_path = os.path.abspath(cargs.output_gfa)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    logger.debug('Subsetting GFA...')
    with open(output_path, 'w') as graph_out:

        with xopen.xopen(cargs.input_gfa, 'r') as graph_in:
            for ln, line in enumerate(graph_in):
                if ln not in keep_lines:
                    continue
                _ = graph_out.write(line)

    logger.debug('GFA saved')

    table_output = os.path.abspath(cargs.output_table)
    os.makedirs(os.path.dirname(table_output), exist_ok=True)
    with open(table_output, 'w') as table:
        _ = table.write('Name,Color,Label\n')
        for tig in select_tigs:
            _, label, color = sorted(annotations[tig], reverse=True)[0]
            _ = table.write('{},{},{}\n'.format(tig, color, label))
    
    logger.debug('CSV table saved')
    return


if __name__ == "__main__":
    logger = None
    rc = 0
    try:
        log_msg_format = "%(asctime)s | %(levelname)s | %(message)s"
        cargs = parse_command_line()
        if cargs.debug:
            log.basicConfig(stream=sys.stdout, level=log.DEBUG, format=log_msg_format)
        else:
            log.basicConfig(stream=sys.stderr, level=log.WARNING, format=log_msg_format)
        logger = log.getLogger()
        logger.debug("Logger initiated")
        main(logger, cargs)
        logger.debug("Run completed - exit")
        log.shutdown()
    except Exception as exc:
        rc = 1
        if logger is not None:
            logger.error("Unrecoverable error: {}".format(str(exc)))
            logger.debug("=== TRACEBACK ===\n\n")
            buf = io.StringIO()
            trb.print_exc(file=buf)
            logger.error(buf.getvalue())
            logger.debug("Exit\n")
            log.shutdown()
        else:
            trb.print_exc()
    finally:
        sys.exit(rc)
