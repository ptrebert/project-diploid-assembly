#!/usr/bin/env python3

__desc__= """

"""

import sys
import os
import argparse
import difflib
import collections as col

import yaml

option_defaults = {
    'repo_folder': os.getcwd(),
    'exec_folder': os.path.split(os.getcwd())[0],
    'sample_name': None,
    'super_population': 'mammals',
    'population': 'humans',
    'family_name': 'unknown',
    'lr_folder': None,
    'lr_project': 'project',
    'lr_seq_platform': None,
    'lr_read_type': None,
    'lr_input_format': None,
    'ss_folder': None,
    'ss_project': 'project',
    'ss_seq_platform': 'ilnxs',
    'ss_read_type': 'npe',
    'ss_library_fractions': None,
}

options_lr_seq_platform = ['pbsq1', 'pbsq2', 'pbrs1', 'pbrs2', 'ontmn', 'ontgd', 'ontpm']
options_lr_read_type = ['clr', 'ccs', 'ul', 'reg']
options_lr_input_format = ['fastq', 'pacbio_native']
options_ss_library_fractions = [1, 2]

def parse_command_line():

    parser = argparse.ArgumentParser(prog='autoconf.py', description=__desc__)

    parser.add_argument(
        '--repository-folder',
        '-rf',
        default=option_defaults['repo_folder'],
        type=str,
        dest='repo_folder',
        help='Full path to folder containing the Snakemake pipeline code. Default: {}'.format(os.getcwd())
    )
    parser.add_argument(
        '--execution-folder',
        '-ef',
        default=option_defaults['exec_folder'],
        type=str,
        dest='exec_folder',
        help='Full path to folder for pipeline execution.'
    )

    sample_group = parser.add_argument_group('Sample configuration')
    sample_group.add_argument(
        '--sample-name',
        '-sn',
        default=option_defaults['sample_name'],
        type=str,
        dest='sample_name',
        help='Name of the sample/individual. It is recommended to use only uppercase letters and numbers.'
    )
    sample_group.add_argument(
        '--super-population',
        '-sp',
        default=option_defaults['super_population'],
        type=str,
        dest='super_population',
        help='Name or ID of super population.'
    )
    sample_group.add_argument(
        '--population',
        '-pop',
        default=option_defaults['population'],
        type=str,
        dest='population',
        help='Name of ID of population.'
    )
    sample_group.add_argument(
        '--family-name',
        '-fam',
        default=option_defaults['family_name'],
        type=str,
        dest='family_name',
        help='Name or ID of family'
    )

    longread_group = parser.add_argument_group('Long-read input configuration')
    longread_group.add_argument(
        '--lr-data-folder',
        '-lrdf',
        default=option_defaults['lr_folder'],
        type=str,
        dest='lr_folder',
        help='Full path to folder containing long-read input data.'
    )
    longread_group.add_argument(
        '--lr-project',
        '-lrp',
        default=option_defaults['lr_project'],
        type=str,
        dest='lr_project',
        help='Project name for long-read input data.'
    )
    longread_group.add_argument(
        '--lr-seq-platform',
        '-lrsp',
        default=option_defaults['lr_seq_platform'],
        type=str,
        choices=options_lr_seq_platform,
        dest='lr_seq_platform',
        help='Specify long-read sequencing platform: {}'.format(options_lr_seq_platform)
    )
    longread_group.add_argument(
        '--lr-read-type',
        '-lrrt',
        default=option_defaults['lr_read_type'],
        type=str,
        choices=options_lr_read_type,
        dest='lr_read_type',
        help='Specify read type: {}'.format(options_lr_read_type)
    )
    longread_group.add_argument(
        '--lr-input-format',
        '-lrif',
        default=option_defaults['lr_input_format'],
        type=str,
        choices=options_lr_input_format,
        help='Specify input data format: {} (FASTQ must be gzipped)'.format(options_lr_input_format)
    )

    strandseq_group = parser.add_argument_group('Strand-seq input configuration')
    strandseq_group.add_argument(
        '--ss-data-folder',
        '-ssdf',
        default=option_defaults['ss_folder'],
        type=str,
        dest='ss_folder',
        help='Full path to folder containing Strand-seq input data (gzipped FASTQ).'
    )
    strandseq_group.add_argument(
        '--ss-project',
        '-ssp',
        default=option_defaults['ss_project'],
        type=str,
        dest='ss_project',
        help='Project name for Strand-seq input data.'
    )
    strandseq_group.add_argument(
        '--ss-seq-platform',
        '-sssp',
        default=option_defaults['ss_seq_platform'],
        type=str,
        dest='ss_seq_platform',
        help='Specify Strand-seq platform.'
    )
    strandseq_group.add_argument(
        '--ss-read-type',
        '-ssrt',
        default=option_defaults['ss_read_type'],
        type=str,
        dest='ss_read_type',
        help='Specify Strand-seq read type.'
    )
    strandseq_group.add_argument(
        '--ss-library-fractions',
        '-sslf',
        default=option_defaults['ss_library_fractions'],
        type=int,
        choices=options_ss_library_fractions,
        dest='ss_library_fractions',
        help='Specify number of different Strand-seq library fractions: {}'.format(options_ss_library_fractions)
    )

def main():


if __name__ == '__main__':
    main()