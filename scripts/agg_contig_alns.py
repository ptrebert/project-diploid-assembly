#!/usr/bin/env python3

import os
import argparse
import functools as fnt
import multiprocessing as mp
import random as rand

import numpy as np
import pandas as pd
import numpy.ma as ma


def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--alignments',
        '-a',
        dest='alignments',
        type=str,
    )
    parser.add_argument(
        '--high-mapq',
        '-hm',
        type=int,
        default=61,
        dest='mapq_high',
        help='Alignment selection is " less / < " than high MAPQ. Default: 61'
    )
    parser.add_argument(
        '--low-mapq',
        '-lm',
        type=int,
        default=0,
        dest='mapq_low',
        help='Alignment selection is " greater or equal / >= " than low MAPQ. Default: 0'
    )
    parser.add_argument(
        '--quantifier',
        '-q',
        choices=['all', 'any'],
        required=True,
        type=str,
        dest='quantifier'
    )
    parser.add_argument(
        '--dump',
        '-d',
        choices=['alignments', 'inverse'],
        type=str,
        default='alignments',
        dest='dump'
    )
    parser.add_argument(
        '--select-chroms',
        '-sc',
        default='"chr[0-9X]+$"',
        type=str,
        dest='select_chroms'
    )
    parser.add_argument(
        '--select-super-pop',
        '-sp',
        nargs='*',
        default=['AFR', 'AMR', 'EAS', 'EUR', 'SAS'],
        dest='select_pop'
    )
    parser.add_argument(
        '--select-tech',
        '-st',
        nargs='*',
        default=['CLR', 'HiFi'],
        dest='select_tech'
    )
    mut_inc_exc = parser.add_mutually_exclusive_group(required=False)
    mut_inc_exc.add_argument(
        '--include-samples',
        '-inc',
        nargs='*',
        default=None,
        dest='include_samples'
    )
    mut_inc_exc.add_argument(
        '--exclude-samples',
        '-exc',
        nargs='*',
        default=None,
        dest='exclude_samples'
    )
    parser.add_argument(
        '--jobs',
        '-j',
        default=1,
        type=int,
        dest='jobs'
    )
    parser.add_argument(
        '--output',
        '-o',
        dest='output',
        type=str
    )
    args = parser.parse_args()
    return args


def keep_alignments(populations, technologies, samples, store_key):
    # store_key = os.path.join(super_pop, sample, technology, haplotype)
    super_pop, sample, tech, _ = store_key.strip('/').split('/')
    return super_pop in populations and tech in technologies and sample in samples


def select_matching_keys(aln_store, select_pop, select_tech, include_samples, exclude_samples):

    aln_keys = []
    
    with pd.HDFStore(aln_store, 'r') as hdf:
        sample_table = hdf[os.path.join('', 'metadata', 'sample')]
        aln_keys = [k for k in hdf.keys() if 'metadata' not in k]

    if 'individual' in sample_table:
        sample_column = 'individual'
    elif 'sample' in sample_table:
        sample_column = 'sample'
    else:
        raise ValueError('Cannot identify sample table in metadata: {}'.format(sample_table.head()))

    sample_table = sample_table.loc[sample_table['2020_SKIP'] != 1, :].copy()

    samples = set(sample_table[sample_column].values)

    if include_samples is None and exclude_samples is None:
        # keep all samples
        select_samples = samples
    elif include_samples is not None:
        select_samples = set(include_samples)
        for s in select_samples:
            if s not in samples:
                raise ValueError('Specified sample does not exist in data set: {}'.format(s))
    elif exclude_samples is not None:
        select_samples = samples - set(exclude_samples)
        if len(select_samples) == 0:
            raise ValueError('No samples left to use after excluding the following: {}'.format(exclude_samples))
    else:
        raise RuntimeError('Unexpected condition')

    keep_key = fnt.partial(keep_alignments, *(select_pop, select_tech, select_samples))
    keys = [k for k in aln_keys if keep_key(k)]

    # limit select_samples to those that match also
    # the other criteria (tech, pop)
    aln_samples = set([ak.strip('/').split('/')[1] for ak in aln_keys])
    select_samples = select_samples.intersect(aln_samples)
    assert len(select_samples) > 0, 'Intersecting samples from aln. store keys with user selection failed'

    return tuple(sorted(keys)), select_samples


def select_matching_chroms(aln_store, chrom_match):

    chrom_match = chrom_match.strip('"')

    with pd.HDFStore(aln_store, 'r') as hdf:
        chroms = hdf[os.path.join('', 'metadata', 'chrom')]
        chroms = chroms.loc[chroms['chrom'].str.match(chrom_match), :].copy()

    return chroms


def build_boolean_mask(chrom, chrom_size, mapq_low, mapq_high, quantifier, aln_store, aln_keys):

    if quantifier == 'all':
        # init: everything TRUE, single case of no
        # alignment to region sufficient to set
        # to FALSE
        mask = np.ones(chrom_size, dtype=np.bool)
        bool_op = np.logical_and
        set_unaligned = True
    elif quantifier == 'any':
        # init: everything FALSE, single case of
        # alignment to region sufficient to set
        # to TRUE
        mask = np.zeros(chrom_size, dtype=np.bool)
        bool_op = np.logical_or
        set_unaligned = False
    else:
        raise ValueError('Unexpected value for quantifier: {}'.format(quantifier))
    
    with pd.HDFStore(aln_store, 'r') as hdf:
        for key in aln_keys:
            alignments = hdf[key]
            select_mapq_high = alignments['mapq'] < mapq_high
            select_mapq_low = alignments['mapq'] >= mapq_low
            select_chrom = alignments['chrom'] == chrom
            alignments = alignments.loc[(select_chrom & select_mapq_low & select_mapq_high), :]
            if alignments.empty:
                continue
            # this stepper is needed to keep track
            # of alignment-free regions per sample;
            # start at beginning of chromosome
            last_end = 0
            # TODO: vec if slow
            for idx, region in alignments.iterrows():
                # if quantifier is ALL, we need to reset regions
                # w/o alignments in this sample that may still
                # be TRUE due to alignments in (all) previous samples
                if set_unaligned and last_end < region['start']:
                    # "last_end < start" needs to be checked for
                    # the case of overlapping alignments
                    mask[last_end:region['start']] = bool_op(mask[last_end:region['start']], False)
                
                mask[region['start']:region['end']] = bool_op(mask[region['start']:region['end']], True)
                # need max here because alignments can be fully contained
                # in other alignments
                last_end = max(last_end, region['end'])
            if set_unaligned:
                mask[last_end:chrom_size] = bool_op(mask[last_end:chrom_size], False)
    
    return mask


def process_alignments(parameter_set):

    chrom, chrom_size, mapq_low, mapq_high, quantifier, dump_type, aln_store, aln_keys = parameter_set

    bool_mask = build_boolean_mask(
        chrom,
        chrom_size,
        mapq_low,
        mapq_high,
        quantifier,
        aln_store,
        aln_keys
    )

    if dump_type == 'alignments':
        # the mask is build by iterating over / masking
        # alignment regions, so need to _invert_
        # the mask to mask out all unaligned
        # (unselected) regions
        genomic_coordinates = ma.masked_array(range(chrom_size), mask=~bool_mask)
    elif dump_type == 'inverse':
        # mask can be used as-is, all unaligned
        # (unselected) regions will be dumped
        genomic_coordinates = ma.masked_array(range(chrom_size), mask=bool_mask)
    else:
        raise ValueError('Unexpected value for "dump_type": {}'.format(dump_type))

    regions = [(s.start, s.stop) for s in  ma.clump_unmasked(genomic_coordinates)]

    return chrom, regions


def build_generic_output_name(dump_regions, mapq_low, mapq_high, quantifier, select_pop, select_tech, select_samples):

    if dump_regions == 'alignments':
        mask = 'MSK:INV'
    elif dump_regions == 'inverse':
        mask = 'MSK:ALN'
    else:
        raise ValueError('Unexpected "dump_regions" parameter: {}'.format(dump_regions))

    mapq = 'LOQ:{}|HIQ:{}'.format(mapq_low, mapq_high)
    quant = 'QNT:{}'.format(quantifier.upper())

    samples = 'SMP:{}'.format(len(select_samples))

    pops = '|'.join(sorted(select_pop))
    techs = '|'.join(sorted(select_tech))
    out_name = '@'.join([mask, mapq, quant, samples, pops, techs])
    return out_name


def main():
    args = parse_command_line()

    alignments_to_process, select_samples = select_matching_keys(
        args.alignments,
        args.select_pop,
        args.select_tech,
        args.include_samples,
        args.exclude_samples
    )

    chroms_to_process = select_matching_chroms(
        args.alignments,
        args.select_chroms
    )

    process_params = []
    for idx, row in chroms_to_process.iterrows():
        param_set = (
            row['chrom'],
            row['size'],
            args.mapq_low,
            args.mapq_high,
            args.quantifier,
            args.dump,
            args.alignments,
            alignments_to_process
        )
        process_params.append(param_set)

    # shuffle once to avoid that all large
    # chromosomes are processed at once
    rand.shuffle(process_params)

    output_regions = []
    out_region_name = build_generic_output_name(
        args.dump,
        args.mapq_low,
        args.mapq_high,
        args.quantifier,
        args.select_pop,
        args.select_tech,
        select_samples
    )
    done = 0
    with mp.Pool(min(args.jobs, len(process_params))) as pool:
        result_iter = pool.imap_unordered(process_alignments, process_params)
        for chrom, regions in result_iter:
            done += 1
            print('{}/{} - finished {}: {} regions'.format(done, len(process_params), chrom, len(regions)))
            subset = pd.DataFrame(regions, columns=['start', 'end'])
            subset['chrom'] = chrom
            subset['name'] = out_region_name
            output_regions.append(subset)

    output_regions = pd.concat(output_regions, axis=0, ignore_index=False)
    output_regions.sort_values(['chrom', 'start', 'end'], inplace=True)

    output_regions = output_regions[['chrom', 'start', 'end', 'name']]
    with open(args.output, 'w') as dump:
        _ = dump.write('#')
        output_regions.to_csv(dump, sep='\t', header=True, index=False)

    return


if __name__ == '__main__':
    rand.seed()
    main()
