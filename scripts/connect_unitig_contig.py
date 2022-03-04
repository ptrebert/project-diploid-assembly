#!/usr/bin/env python3

from multiprocessing.sharedctypes import Value
import sys
import argparse as argp
import traceback as trb
import pathlib as pl
import hashlib
import collections as col

import pandas as pd
import networkx as nx


def parse_command_line():
    parser = argp.ArgumentParser()
    parser.add_argument(
        "--unitig-gfa",
        "-u",
        type=lambda x: pl.Path(x).resolve().absolute(),
        dest="unitigs"
    )
    parser.add_argument(
        "--contig-gfa",
        "-c",
        type=lambda x: pl.Path(x).resolve().absolute(),
        dest="contigs"
    )
    parser.add_argument(
        "--saarclusters",
        "-s",
        type=lambda x: pl.Path(x).resolve().absolute(),
        dest="saarclusters"
    )
    parser.add_argument(
        "--min-contig-read-depth",
        "-rd",
        type=int,
        default=5,
        dest="contig_read_depth"
    )
    parser.add_argument(
        "--ignore-non-assembly-sequence",
        "-inas",
        action="store_true",
        default=False,
        dest="ignore_extra_seq",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=lambda x: pl.Path(x).resolve().absolute(),
        dest="output"
    )
    parser.add_argument(
        "--unassigned-unitigs",
        "-uu",
        type=lambda x: pl.Path(x).resolve().absolute(),
        dest="out_unassigned",
        default=None
    )
    args = parser.parse_args()
    return args


def process_graph_file(file_path, keep_edges):

    tig_to_readcount = col.Counter()
    tig_to_length = col.Counter()
    tig_to_read_depth = col.Counter()
    read_to_tig = dict()
    edges = []
    nodes = []

    keep_records = ['S', 'A']
    if keep_edges:
        keep_records.append('L')

    with open(file_path, 'r') as graph:
        for ln, line in enumerate(graph, start=1):
            if not line[0] in keep_records:
                continue
            parts = line.strip().split()
            if parts[0] == 'S':
                tig = parts[1]
                tig_length = int(parts[-2].split(':')[-1])
                tig_read_depth = int(parts[-1].split(':')[-1])
                tig_to_length[tig] = tig_length
                tig_to_read_depth[tig] = tig_read_depth
                nodes.append(tig)
            elif parts[0] == 'A':
                tig = parts[1]
                read = parts[4]
                if read in read_to_tig:
                    seen_tig = read_to_tig[read]
                    # can happen that read is used in frw. and rev.
                    assert seen_tig == tig
                read_to_tig[read] = tig
                tig_to_readcount[tig] += 1
            elif parts[0] == 'L':
                tig1 = parts[1]
                tig2 = parts[3]
                edges.append((tig1, tig2))
            else:
                RuntimeError(f'GFA parser error: {ln} / {line.strip()}')
    return tig_to_readcount, tig_to_length, tig_to_read_depth, read_to_tig, edges, nodes


def process_raw_unitig_graph(file_path):

    tig_to_readcount, tig_to_length, tig_to_read_depth, read_to_tig, edges, nodes = process_graph_file(file_path, True)
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    graph.add_edges_from(edges)

    tig_to_cc_hash = dict()
    tig_to_cc_size = dict()
    for cc_nodes in nx.connected_components(graph):
        cc_id = ''.join(sorted(cc_nodes))
        cc_key = hashlib.md5(cc_id.encode('ascii')).hexdigest()
        for node in cc_nodes:
            tig_to_cc_hash[node] = cc_key
            tig_to_cc_size[node] = len(cc_nodes)

    unitig_table = pd.DataFrame.from_records(
        [tig_to_readcount, tig_to_cc_size, tig_to_cc_hash, tig_to_length, tig_to_read_depth],
        index=['unitig_read_count', 'component_size', 'component_hash', 'unitig_length', 'unitig_read_depth']
    )
    unitig_table = unitig_table.transpose()
    return unitig_table, read_to_tig


def build_unitig_table(args):
    """
    Potential source of confusion:
    In rare cases (e.g., utg022947l and utg025333l for HG02257), a unitig can be
    assigned to several primary contigs on the basis of the read distribution between
    unitig and contig. Since a 1-to-1 mapping seems desirable for the implemented
    heuristic of patching back in assembled sequence that is discarded by SaaRclust,
    the below procedure determines the assignment following a simple majority vote.

    2022-03-04 --- sample HG02257
    The above situation now triggers an exception because the contig ptg000250l is
    not assigned to any unitig (reads m64076_200129_001835/85721750/ccs and
    m64076_200130_064345/127666640/ccs establish relation, but more reads of the
    unitig utg022947l relate to another contig [ptg000188l]). Hence, the contig
    ptg000250l is detected as an unknown sequence, which is misleading.

    """
    unitig_table, read_to_unitig = process_raw_unitig_graph(args.unitigs)
    # skip over contig length - part of SaaRclust output
    contig_readcount, _, contig_read_depth, read_to_contig, _, _ = process_graph_file(args.contigs, False)

    unitigs_to_contig_count = col.defaultdict(col.Counter)
    contigs_to_unitigs_count = col.defaultdict(set)
    unitig_to_contig = dict()
    unseen_reads = col.Counter()
    for read, contig in read_to_contig.items():
        try:
            unitig = read_to_unitig[read]
        except KeyError:
            unseen_reads[contig] += 1
        # count unitig to contig assignments - there is multiplicity
        # because we are iterating over reads
        unitigs_to_contig_count[unitig][contig] += 1
        contigs_to_unitigs_count[contig].add(unitig)
    contigs_to_unitigs_count = col.Counter({c: len(u) for c, u in contigs_to_unitigs_count.items()})

    # make majority decision explicit here
    unitig_to_contig = dict()
    minority_contigs = set()
    majority_contigs = set()
    for unitig, contig_counts in unitigs_to_contig_count.items():
        contig, _ = contig_counts.most_common(1)[0]
        # fix here [2022-03-04]
        # keep a record of minority vote contigs,
        # if any left (= unassigned) at the end,
        # fix that in the table below
        majority_contigs.add(contig)
        [minority_contigs.add(contig) for contig, _ in contig_counts.most_common()[1:]]
        unitig_to_contig[unitig] = contig

    unassigned_contigs = minority_contigs - majority_contigs

    unitig_table['unitig'] = unitig_table.index
    unitig_table.reset_index(drop=True, inplace=True)
    unitig_table['contig'] = unitig_table['unitig'].map(lambda x: unitig_to_contig.get(x, 'unassigned'))
    if unassigned_contigs:
        records = []
        for ctg in unassigned_contigs:
            record = {
                'contig': ctg,
                'unitig': 'unassigned',
                'component_hash': 'unknown',
                'unitig_read_count': 0,
                'component_size': 0,
                'unitig_length': 0,
                'unitig_read_depth': 0
            }
            records.append(record)
        tmp = pd.DataFrame.from_records(records)
        unitig_table = pd.concat([unitig_table, tmp], axis=0, ignore_index=False)

    unitig_table['contig_read_count'] = unitig_table['contig'].map(lambda x: contig_readcount[x])
    unitig_table['contig_read_depth'] = unitig_table['contig'].map(lambda x: contig_read_depth[x])
    unitig_table['contig_only_reads'] = unitig_table['contig'].map(lambda x: unseen_reads[x])
    unitig_table['contig_unitig_count'] = unitig_table['contig'].map(lambda x: contigs_to_unitigs_count[x])
    unitig_table['unitig'] = unitig_table.index
    unassigned_unitigs = unitig_table.loc[unitig_table['contig'] == 'unassigned', :].copy()
    unitig_table = unitig_table.loc[~unitig_table.index.isin(unassigned_unitigs.index), :].copy()
    unitig_table.reset_index(drop=True, inplace=True)
    unassigned_unitigs.reset_index(drop=True, inplace=True)
    return unitig_table, unassigned_unitigs


ORIENT_MAP = {
    'dir': 'F',
    'revcomp': 'R'
}


def derive_scl_key(row):

    selected = 'sY' if row['is_selected'] else 'sN'
    clustered = 'cY' if row['is_clustered'] else 'cN'
    error = 'eY' if row['has_error'] else 'eN'
    direction = 'd' + ORIENT_MAP.get(row['has_direction'], 'U')
    ploidy = 'p'
    if row['ploidy_estimate'] in ['1n', '2n']:
        ploidy += row['ploidy_estimate'].strip('n')
    elif pd.isna(row['ploidy_estimate']):
        ploidy += '0'
    else:
        raise ValueError(f'Cannot process ploidy in SaaRclust row: {row}')
    scl_key = f'{selected}{clustered}{error}{direction}{ploidy}'
    return scl_key


def make_num_cluster_id(cluster_id):
    if pd.isna(cluster_id):
        cid = 0
    else:
        cid = int(cluster_id.strip('cluster'))
    return cid


def make_padded_cluster_id(cluster_id):
    if pd.isna(cluster_id):
        pid = '00'
    else:
        cid = cluster_id.strip('cluster')
        pid = f'{cid:0>2}'
    return pid


def process_saarclust_table(file_path):

    default_header = [
        "ctg",
        "ctg.len",
        "sizeSelect",
        "Clustered",
        "putative.error",
        "Dir",
        "Cluster.ID",
        "Ploidy"
    ]
    saarclusters = pd.read_csv(file_path, sep='\t', header=0)
    if saarclusters.columns.tolist() != default_header:
        raise RuntimeError(f'Format change in SaaRclust output table: {saarclusters.columns}')
    new_header = [
        'contig', 'contig_length', 'is_selected', 'is_clustered',
        'has_error', 'has_direction', 'cluster_id', 'ploidy_estimate'
    ]
    saarclusters.columns = new_header
    saarclusters['scl_key'] = saarclusters.apply(derive_scl_key, axis=1)
    saarclusters['has_direction'].fillna('unk', inplace=True)
    saarclusters['cluster_id'].fillna('cluster99', inplace=True)
    saarclusters['ploidy_estimate'].fillna('0n', inplace=True)
    return saarclusters


def run_guilt_by_association(unitig_table, contig_read_depth):

    ignore_low_read_depth = unitig_table['contig_read_depth'] >= contig_read_depth
    process_unassigned = unitig_table['cluster_id'] == 'cluster99'
    select_unassigned_contigs = ignore_low_read_depth & process_unassigned

    reassign = []

    unassigned_contigs = unitig_table.loc[select_unassigned_contigs, 'contig'].unique()
    for ctg in unassigned_contigs:
        # check neighborhood
        components = unitig_table.loc[unitig_table['contig'] == ctg, 'component_hash'].unique().tolist()
        selector = (unitig_table['component_hash'].isin(components)) & (unitig_table['cluster_id'] != 'cluster99')
        subgraph = unitig_table.loc[selector, :].drop_duplicates(['contig', 'cluster_id'])
        subgraph_clusters = subgraph.groupby('cluster_id')['contig_length'].sum()
        if subgraph_clusters.shape[0] != 1:
            continue
        reassign.append((ctg, subgraph_clusters.index[0]))
    unitig_table['init_cluster_id'] = unitig_table['cluster_id']
    if reassign:    
        for ctg, new_cluster in reassign:
            unitig_table.loc[unitig_table['contig'] == ctg, 'cluster_id'] = new_cluster
    return unitig_table

def reorder_table(unitig_table):

    header = [
        'cluster_id',
        'contig',
        'contig_length',
        'contig_read_depth',
        'contig_read_count',
        'contig_only_reads',
        'contig_unitig_count',
        'unitig',
        'unitig_length',
        'unitig_read_count',
        'unitig_read_depth',
        'component_size',
        'component_hash',
        'scl_key',
        'cluster_num_id',
        'cluster_pad_id',
        'init_cluster_id',
        'is_selected',
        'is_clustered',
        'has_error',
        'has_direction',
        'ploidy_estimate'
    ]
    assert all(c in header for c in unitig_table.columns)
    unitig_table = unitig_table[header]
    # these columns included NaN (float) at some point
    unitig_table['component_size'] = unitig_table['component_size'].astype(int)
    unitig_table['unitig_length'] = unitig_table['unitig_length'].astype(int)
    unitig_table['unitig_read_count'] = unitig_table['unitig_read_count'].astype(int)
    unitig_table['unitig_read_depth'] = unitig_table['unitig_read_depth'].astype(int)
    return unitig_table[header]


def main():
    args = parse_command_line()
    unitig_table, unassigned_unitigs = build_unitig_table(args)
    saarclusters = process_saarclust_table(args.saarclusters)
    unitig_table = unitig_table.merge(saarclusters, on='contig', how='outer')
    # change here: if parts (contigs) of the assembly are replaced by new contigs,
    # the original contig names will not show up in SaaRclust's contig report, and are hence
    # unknown in this context. Since this script reads to original
    # hifiasm graph, the original contig names are reintroduced here, causing a mismatch
    unclustered_contigs = pd.isnull(unitig_table['is_selected'])
    if not unitig_table.loc[unclustered_contigs, :].empty:
        if args.ignore_extra_seq:
            # drop the reintorduced contig names from the table
            unitig_table = unitig_table.loc[~unclustered_contigs, :].copy()
        else:
            raise ValueError(f'Unknown sequence detected in NHR assembly:\n=====\n{unitig_table.loc[unclustered_contigs, :]}\n=====\n')
    # change here: unassigned unitigs were discarded before because of the default "inner"
    # merge; does not affect downstream process, but to have the full assembly record
    # available, keep them in a separate output file
    no_unitig = pd.isnull(unitig_table['unitig'])
    yes_contig = ~pd.isnull(unitig_table['contig'])
    unknown_sequence = no_unitig & yes_contig
    if not unitig_table.loc[unknown_sequence, :].empty:
        if args.ignore_extra_seq:
            int_columns = [
                'unitig_read_count',
                'component_size',
                'contig_length',
                'unitig_length',
                'unitig_read_depth',
                'contig_read_count',
                'contig_read_depth',
                'contig_only_reads',
                'contig_unitig_count'
            ]
            str_columns = ['component_hash', 'unitig']
            for int_column in int_columns:
                unitig_table[int_column].fillna(0, inplace=True)
                unitig_table[int_column] = unitig_table[int_column].astype(int)

            for str_column in str_columns:
                unitig_table[str_column] = 'unknown'
        else:
            raise ValueError(f'Unknown sequence detected in clustered assembly:\n=====\n{unitig_table.loc[unknown_sequence, :]}\n=====\n')
    unitig_table = run_guilt_by_association(unitig_table, args.contig_read_depth)
    unitig_table['cluster_num_id'] = unitig_table['cluster_id'].map(make_num_cluster_id)
    unitig_table.sort_values(['cluster_num_id', 'contig_length'], ascending=[True, False], inplace=True)
    unitig_table['cluster_pad_id'] = unitig_table['cluster_id'].map(make_padded_cluster_id)

    assert not pd.isna(unitig_table).any(axis=1).any()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    unitig_table = reorder_table(unitig_table)
    unitig_table.to_csv(args.output, sep='\t', header=True, index=False)

    if args.out_unassigned is not None:
        args.out_unassigned.parent.mkdir(parents=True, exist_ok=True)
        unassigned_unitigs.to_csv(args.out_unassigned, sep='\t', header=True, index=False)

    return 0


if __name__ == '__main__':
    sys.exit(main())
