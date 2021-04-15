#!/usr/bin/env python3

import os as os
import sys as sys
import io as io
import traceback as trb
import argparse as argp
import logging as log

# change in the numexpr package / dep of pandas
# if not set, issues pointless info logging message
os.environ['NUMEXPR_MAX_THREADS'] = '2'

import numpy as np
import pandas as pd


__docs__ = """
Script to read GAF alignment files produced by GraphAligner.
Script developed for ONT to HiFi graph alignment.
Selects single best alignment per ONT read (max alignment score).
"""

GAF_HEADER = [
    'read_name',
    'read_length',
    'read_align_start',
    'read_align_end',
    'read_align_orientation',
    'path',
    'path_length',
    'path_align_start',
    'path_align_end',
    'num_matches',
    'align_block_length',
    'mapq',
    'edit_distance',
    'alignment_score',
    'divergence',
    'identity'
]

GAF_COLUMN_TYPES = {
    'read_name': str,
    'read_length': np.uint32,
    'read_align_start': np.uint32,
    'read_align_end': np.uint32,
    'read_align_orientation': str,
    'path': str,
    'path_length': np.uint32,
    'path_align_start': np.uint32,
    'path_align_end': np.uint32,
    'num_matches': np.uint32,
    'align_block_length': np.uint32,
    'mapq': np.uint8,
    #'edit_distance': np.uint32,  -- see converter
    #'alignment_score': np.uint32,  -- see converter
    'divergence': str,
    #'identity': np.uint8  -- see converter
}

IGNORE_COLUMNS = [
    'read_align_orientation',
    'path_align_start',
    'path_align_end',
    'num_matches',
    'align_block_length',
    'mapq',
    'divergence'
]

USE_COLUMNS = [c for c in GAF_HEADER if c not in IGNORE_COLUMNS]


def _convert_nm_column(nm_tag):
    """
    Convert NM tag / edit distance
    NM:i:3269
    """
    return np.uint32(nm_tag.split(':')[-1])


def _convert_as_column(as_tag):
    """
    Convert AS tag / alignment score
    Round to int
    AS:f:24717.1
    """
    return np.uint32(round(float(as_tag.split(':')[-1]), 0))
    

def _convert_id_column(id_tag):
    """
    Convert id tag / pct. sequence identity
    Round to int to save some space
    (1 - dv) tag
    """
    return np.uint8(round(float(id_tag.split(':')[-1]) * 100, 0))


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog="filter_gaf.py", description=__doc__)
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
        "--input",
        "-i",
        type=str,
        required=True,
        dest="input",
        help="Full path GAF (can be compressed).",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        dest="output",
        help="Full path to output file (pandas/pytables HDF format). Directories will "
             "be created if they do not exist. File name must have extension '.h5'!"
    )
    parser.add_argument(
        "--size-fractions",
        "-sf",
        type=int,
        default=None,
        dest="fractions",
        nargs="*",
        help="Dump read and path information for these read length fractions. "
            "Specify list of lower boundaries, greater-or-equal will be dumped."
    )
    return parser.parse_args()

def read_gaf_file(input_file, logger):
    """
    """
    logger.debug('Reading GAF file: {}'.format(input_file))
    df = pd.read_csv(
        input_file,
        sep='\t',
        header=None,
        index_col=False,
        names=GAF_HEADER,
        usecols=USE_COLUMNS,
        dtype=GAF_COLUMN_TYPES,
        converters={
            'edit_distance': _convert_nm_column,
            'alignment_score': _convert_as_column,
            'identity': _convert_id_column
        },
        encoding='ascii',
        low_memory=True
    )
    logger.debug('Reading complete')

    df['edit_distance'] = df['edit_distance'].astype(np.uint32)
    df['alignment_score'] = df['alignment_score'].astype(np.uint32)
    df['identity'] = df['identity'].astype(np.uint8)

    df['read_align_length'] = df['read_align_end'] - df['read_align_start']
    df['read_align_length'] = df['read_align_length'].astype(np.uint32)

    total_alignments = df.shape[0]
    logger.debug('Raw alignments in GAF: {}'.format(total_alignments))

    df.set_index(np.arange(1, df.shape[0]+1, dtype=np.uint32), inplace=True)
    df.index.name = 'gaf_line_number'  # important if the original GAF file should also be reduced later

    logger.debug('Dropping duplicate read IDs...')
    
    # if selection strategy changes / not on the basis of max alignment score,
    # can enumerate duplicate read IDs like this
    #df['enum'] = df.groupby('read_name').cumcount()

    # keep first: relies on GraphAligner behavior to write best alignment first
    df.drop_duplicates('read_name', keep='first', ignore_index=False, inplace=True)

    uniq_alignments = df.shape[0]
    logger.debug('GAF alignments w/o duplicates: {}'.format(uniq_alignments))
    multi_aln_ratio = round(total_alignments / uniq_alignments, 1)
    logger.debug('Mean alignments per read: {}'.format(multi_aln_ratio))

    return df


def store_gaf_as_hdf(gaf, output_file, logger):
    """
    """
    full_output_path = os.path.abspath(output_file)
    out_dir = os.path.dirname(full_output_path)
    logger.debug('Creating output directory: {}'.format(out_dir))
    os.makedirs(out_dir, exist_ok=True)

    gaf.to_hdf(
        output_file,
        'gaf',
        complevel=9,
        encoding='ascii',
        format='fixed'
    )
    return full_output_path


def split_paths(paths):
    """
    """
    translation_table = dict((i,i) for i in '1234567890')
    translation_table['>'] = ' '
    translation_table['<'] = ' '
    translation_table = str.maketrans(translation_table)

    all_paths = paths.apply(lambda x: x.translate(translation_table).split()).values

    nodes = sorted(
        set().union(*all_paths)
    )

    return nodes


def dump_size_fraction_info(gaf, threshold, output_path, logger):
    """
    """
    base_out, h5_ext = output_path.rsplit('.', 1)
    assert h5_ext == 'h5', 'Output file extension should be .h5'

    ec_name_infos = ['read_name', 'read_align_start', 'read_align_end']
    if 'readec_name' not in gaf:
        logger.debug('Setting read[ec] name...')
        # note hard-coded zero because of GraphAligner's alignment output sort-order
        gaf['readec_name'] = gaf[ec_name_infos].apply(lambda x: '{}_0_{}_{}'.format(*x), axis=1)
    
    sub = gaf.loc[gaf['read_align_length'] >= threshold, :]
    logger.debug('Selected {} reads to dump as size fraction greater-or-equal to {}'.format(sub.shape[0], threshold))

    name_infix = 'geq{}'.format(threshold)
    
    logger.debug('Dumping read[ec] names...')
    readec_name_file = '.'.join([base_out, name_infix, 'read-ec', 'txt'])
    sub.to_csv(
        readec_name_file,
        sep='\t',
        header=False,
        index=False,
        encoding='ascii',
        columns=['readec_name']
    )

    logger.debug('Dumping path node IDs...')
    path_node_file = '.'.join([base_out, name_infix, 'path-nodes', 'txt'])
    path_nodes = split_paths(sub['path'])
    with open(path_node_file, 'w', encoding='ascii') as dump:
        _ = dump.write('\n'.join(path_nodes) + '\n')
    
    logger.debug('Dumping read[ec] name to path mapping...')
    readec_path_file = '.'.join([base_out, name_infix, 'readec-path', 'tsv'])
    sub.to_csv(
        readec_path_file,
        sep='\t',
        header=False,
        index=False,
        encoding='ascii',
        columns=['readec_name', 'path']
    )
    logger.debug('Dumped all size fraction infos')
    return



def main(logger, cargs):
    """
    :param logger:
    :param cargs:
    :return:
    """
    if not cargs.output.endswith('.h5'):
        raise ValueError('Please set the output file name to [<SOME-PATH>/]<SOME-NAME>.h5')

    logger.debug('Starting GAF filtering...')

    gaf = read_gaf_file(cargs.input, logger)

    logger.debug('Storing GAF information as HDF')
    output_file = store_gaf_as_hdf(gaf, cargs.output, logger)
    logger.debug('GAF information saved')

    dump_size_fraction_info(
        gaf,
        0,
        output_file,
        logger
    )

    if cargs.fractions is not None:
        logger.debug('Dumping read and path info for size fractions')
        for threshold in cargs.fractions:
                if threshold == 0:
                    continue
                dump_size_fraction_info(
                    gaf,
                    threshold,
                    output_file,
                    logger
                )

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
