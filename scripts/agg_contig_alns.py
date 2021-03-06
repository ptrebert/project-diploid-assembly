#!/usr/bin/env python3

import os
import sys
import argparse
import functools as fnt
import collections as col
import multiprocessing as mp
import random as rand
import math

import numpy as np
import pandas as pd
import numpy.ma as ma


def parse_command_line():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--alignments',
        '-a',
        dest='alignments',
        type=str,
    )
    parser.add_argument(
        '--high-mapq',
        '-hiq',
        type=int,
        default=61,
        dest='mapq_high',
        help='Alignment selection is [less / <] than high MAPQ. Default: 61'
    )
    parser.add_argument(
        '--low-mapq',
        '-loq',
        type=int,
        default=0,
        dest='mapq_low',
        help='Alignment selection is [greater or equal / >=] than low MAPQ. Default: 0'
    )
    parser.add_argument(
        '--quantifier',
        '-q',
        choices=['all', 'any', 'lt', 'geq'],
        required=True,
        type=str,
        dest='quantifier',
        help='Possible selection strategies for alignments between "--low-mapq" and "--high-mapq":\n'
            '*all*: only mark regions where all selected assemblies have alignments.\n'
            '*any*: only mark regions where any of the selected assemblies has an alignment.\n'
            '*lt*: only mark regions where [less / <] than "--select-assemblies" number/percent\n'
            'of the selected assemblies have alignments.\n'
            '*geq*: only mark regions where [greater or equal / >=] than "--select-assemblies" number/percent\n'
            'of the selected assemblies have alignments.'
    )
    parser.add_argument(
        '--select-assemblies',
        '-sa',
        type=float,
        default=-1,
        dest='select_assemblies',
        help='Specify percentage or number of assemblies to threshold for "--quantifier":\n'
            'The default value of "-1" selects all assemblies after modifying the set with\n'
            'samples to be included or excluded, selecting for sex, technology etc (see below).\n'
            'For "--select-value=absolute", the following conditions are checked:\n'
            '*lt*: select region if 0 < num_alignments < num_assemblies\n'
            '*geq*: select region if num_alignments >= num_assemblies\n'
            'For "--select-value=fraction", the following conditions are checked:\n'
            '*lt*: select region if 0 < num_alignments < ceiling(assembly-fraction * num_assemblies)\n'
            '*geq*: select region if num_alignments >= floor(assembly-fraction * num_assemblies)'
    )
    parser.add_argument(
        '--select-value',
        '-sv',
        type=str,
        default='absolute',
        choices=['fraction', 'absolute'],
        dest='select_value',
        help='Specify how the value for "--select-assemblies" should be interpreted.'
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
        '--dump-counts',
        '-dc',
        action='store_true',
        default=False,
        dest='dump_counts',
        help='Dump counts / quantitative mask per chromosome (uses output file name).'
    )
    parser.add_argument(
        '--min-size',
        '-ms',
        default=0,
        type=int,
        dest='min_size',
        help='Minimum region size to include in output [greater or equal / >=].'
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
    parser.add_argument(
        '--tech-priority',
        '-tp',
        default=None,
        choices=[None, 'CLR', 'HiFi'],
        dest='tech_priority',
        help='For samples with alignments from two different technologies, select only the tech with higher priority.'
    )
    parser.add_argument(
        '--select-sex',
        '-sx',
        nargs='*',
        default=['male', 'female'],
        dest='select_sex'
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


def keep_alignments(populations, sexes, technologies, samples, store_key):
    # store_key = os.path.join(super_pop, sample, technology, haplotype)
    super_pop, sample, sex, tech, _ = store_key.strip('/').split('/')

    keep_pop = super_pop in populations
    keep_sex = sex in sexes
    keep_sample = sample in samples
    keep_tech = tech in technologies

    return  keep_pop & keep_sex & keep_sample & keep_tech


def filter_low_priorty_tech(aln_keys, tech_priority):

    sample_counts = col.Counter([k.split('/')[2] for k in aln_keys])
    tech_dups = set()
    for n, c in sample_counts.most_common():
        assert c == 2 or c == 4, 'Unexpected sample count: {} / {}'.format(n, c)
        if c == 2:
            # two haplotypes
            break
        tech_dups.add(n)
    
    filtered_keys = []
    if tech_dups:
        for k in aln_keys:
            super_pop, sample, sex, tech, _ = k.strip('/').split('/')
            if sample in tech_dups and tech != tech_priority:
                continue
            filtered_keys.append(k)
    else:
        filtered_keys = aln_keys
    return filtered_keys


def select_matching_keys(aln_store, select_pop, select_sex, select_tech, tech_priority, include_samples, exclude_samples):

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

    keep_key = fnt.partial(keep_alignments, *(select_pop, select_sex, select_tech, select_samples))
    keys = [k for k in aln_keys if keep_key(k)]
 
    if tech_priority is not None:
        keys = filter_low_priorty_tech(keys, tech_priority)

    aln_samples = set([k.strip('/').split('/')[1] for k in keys])

    # limit select_samples to those that match also
    # the other criteria (tech, pop)
    select_samples = select_samples.intersection(aln_samples)

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


def build_quantitative_mask(chrom, chrom_size, mapq_low, mapq_high, quantifier, assembly_threshold, dump_counts, aln_store, aln_keys):

    quant_mask = np.zeros(chrom_size, dtype=np.int8)
    mask = np.zeros(chrom_size, dtype=np.bool)
    
    with pd.HDFStore(aln_store, 'r') as hdf:
        for key in aln_keys:
            mask &= False
            alignments = hdf[key]
            select_mapq_high = alignments['mapq'] < mapq_high
            select_mapq_low = alignments['mapq'] >= mapq_low
            select_chrom = alignments['chrom'] == chrom
            alignments = alignments.loc[(select_chrom & select_mapq_low & select_mapq_high), :]
            if alignments.empty:
                continue

            for idx, region in alignments.iterrows():
                mask[region['start']:region['end']] = True
            
            quant_mask[mask] += 1
    
    if dump_counts is not None:
        qmask_dump = dump_counts + '.qmask.{}.npy'.format(chrom)
        np.save(qmask_dump, quant_mask, allow_pickle=False)

    mask &= False
    if quantifier == 'lt':
        mask[np.logical_and(0 < quant_mask, quant_mask < assembly_threshold)] = True
    elif quantifier == 'geq':
        mask[(quant_mask >= assembly_threshold)] = True
    else:
        raise ValueError('Unexpected quantifier: {}'.format(quantifier))
    return mask


def process_alignments(parameter_set):

    chrom, chrom_size, mapq_low, mapq_high, quantifier, assembly_threshold, dump_type, dump_counts, aln_store, aln_keys = parameter_set

    if quantifier in ['all', 'any']:
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
            genomic_coordinates = ma.masked_array(np.arange(chrom_size, dtype=np.int32), mask=~bool_mask)
        elif dump_type == 'inverse':
            # mask can be used as-is, all unaligned
            # (unselected) regions will be dumped
            genomic_coordinates = ma.masked_array(np.arange(chrom_size, dtype=np.int32), mask=bool_mask)
        else:
            raise ValueError('Unexpected value for "dump_type": {}'.format(dump_type))
    elif quantifier in ['lt', 'geq']:
        quant_mask = build_quantitative_mask(
            chrom,
            chrom_size,
            mapq_low,
            mapq_high,
            quantifier,
            assembly_threshold,
            dump_counts,
            aln_store,
            aln_keys
        )
        # negating mask: see reason above
        genomic_coordinates = ma.masked_array(np.arange(chrom_size, dtype=np.int32), mask=~quant_mask)
    else:
        raise ValueError('Unsupported quantifier: {}'.format(quantifier))

    regions = [(s.start, s.stop) for s in  ma.clump_unmasked(genomic_coordinates)]

    return chrom, regions


def build_generic_output_name(dump_regions, mapq_low, mapq_high, quantifier, select_pop, select_sex, select_tech, select_samples, assembly_threshold):

    if dump_regions == 'alignments':
        mask = 'MSK:INV'
    elif dump_regions == 'inverse':
        mask = 'MSK:ALN'
    else:
        raise ValueError('Unexpected "dump_regions" parameter: {}'.format(dump_regions))

    mapq = 'LOQ:{}|HIQ:{}'.format(mapq_low, mapq_high)
    quant = 'QNT:{}'.format(quantifier.upper())

    if quantifier in ['lt', 'geq']:
        quant += ':{}'.format(assembly_threshold)

    samples = 'SMP:{}'.format(len(select_samples))

    pops = '|'.join(sorted(select_pop))
    techs = '|'.join(sorted(select_tech))
    sexs = '|'.join(sorted(select_sex))
    out_name = '@'.join([mask, mapq, quant, samples, pops, sexs, techs])
    return out_name


def compute_assembly_threshold(alignments_to_process, quantifier, select_assemblies, select_value):

    if select_assemblies == -1:
        # select all assembly alignments that are left
        # after modifying the set with include/exclude
        # (has happened above)
        assembly_threshold = alignments_to_process
    elif select_assemblies == 0:
        raise ValueError('Parameter "select assemblies" is set to 0. '
                        'If you want select unaligned regions, use the quantifier "any" '
                        'and dump the "inverse" set of marked regions.')
    elif 0 < select_assemblies < 1 and select_value == 'absolute':
        raise ValueError('Parameter "select assemblies" is set to {} but should be '
                        'interpreted as absolute value. Set "select assemblies" to '
                        'at least 1 or set "--select-value" to "fraction".'.format(select_assemblies))
    elif 1 <= select_assemblies and select_value == 'absolute':
        assembly_threshold = int(select_assemblies)
        if alignments_to_process < assembly_threshold:
            raise ValueError('Fewer alignments to process ({}) than selected assemblies ({}).'.format(alignments_to_process, assembly_threshold))
    elif 0 < select_assemblies <= 1 and select_value == 'fraction':
        if quantifier in ['lt', 'geq']:
            if quantifier == 'lt':
                assembly_threshold = math.ceil(alignments_to_process * select_assemblies)
            else:
                assembly_threshold = math.floor(alignments_to_process * select_assemblies)
            assert assembly_threshold < 255, 'Too many assemblies for thresholding'
    else:
        assembly_threshold = int(select_assemblies)

    assert 0 < assembly_threshold, \
        'Computing assembly threshold failed: {} / {} / {} / {}'.format(assembly_threshold, alignments_to_process, select_assemblies, select_value)

    return assembly_threshold

def main():
    args = parse_command_line()

    alignments_to_process, select_samples = select_matching_keys(
        args.alignments,
        args.select_pop,
        args.select_sex,
        args.select_tech,
        args.tech_priority,
        args.include_samples,
        args.exclude_samples
    )

    chroms_to_process = select_matching_chroms(
        args.alignments,
        args.select_chroms
    )

    assembly_threshold = compute_assembly_threshold(
        len(alignments_to_process),
        args.quantifier,
        args.select_assemblies,
        args.select_value
    )

    process_params = []
    for idx, row in chroms_to_process.iterrows():
        if args.dump_counts:
            dump_counts = args.output.rsplit('.', 1)[0]
        else:
            dump_counts = None
        param_set = (
            row['chrom'],
            row['size'],
            args.mapq_low,
            args.mapq_high,
            args.quantifier,
            assembly_threshold,
            args.dump,
            dump_counts,
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
        args.select_sex,
        args.select_tech,
        select_samples,
        assembly_threshold
    )
    done = 0
    with mp.Pool(min(args.jobs, len(process_params))) as pool:
        result_iter = pool.imap_unordered(process_alignments, process_params)
        for chrom, regions in result_iter:
            done += 1
            print('{}/{} - finished {}: {} regions'.format(done, len(process_params), chrom, len(regions)))
            subset = pd.DataFrame(regions, columns=['start', 'end'])
            region_size_select = (subset['end'] - subset['start']) >= args.min_size
            if not region_size_select.any():
                print('No regions left after size thresholding: {}'.format(args.min_size))
                continue
            subset = subset.loc[region_size_select, :].copy()
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
