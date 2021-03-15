#!/usr/bin/env python3

__desc__ = """
autoconf.py creates a configuration YAML file for a single human sample pipeline run
with default parameters. The autoconf.py script assumes a simple run with no
precomputed files being "injected" into the pipeline via file linking.
"""

import os
import argparse
import difflib
import collections as col
import itertools
import hashlib
import warnings

import yaml

option_defaults = col.OrderedDict([
    ('repo_folder', os.getcwd()),
    ('exec_folder', os.path.split(os.getcwd())[0]),
    ('sample_name', None),
    ('super_population', 'animals'),
    ('population', 'mammals'),
    ('family_name', 'humans'),
    ('lr_folder', None),
    ('lr_project', 'project'),
    ('lr_seq_platform', None),
    ('lr_read_type', None),
    ('lr_input_format', None),
    ('ss_folder', None),
    ('ss_project', 'project'),
    ('ss_seq_platform', 'ilany'),
    ('ss_read_type', 'npe'),
    ('ss_library_fractions', None),
    ('param_ver', 13),
    ('max_cpu', None),
    ('high_cpu', None),
    ('medium_cpu', None),
    ('low_cpu', None),
    ('local_copy', 'no'),
    ('singularity_module', ''),
    ('nhr_assembler', None),
    ('hap_assembler', None),
    ('var_caller', None)
])

option_constraints = {
    'lr_seq_platform': {
        'check': ['pbsq1', 'pbsq2', 'pbrs1', 'pbrs2', 'ontmn', 'ontgd', 'ontpm'],
        'display': ['(PacBio Sequel/RS) >> ', 'pbsq1', 'pbsq2', 'pbrs1', 'pbrs2',
                    '(ONT Min/Grid/Prometh-ION) >>', 'ontmn', 'ontgd', 'ontpm'],
    },
    'lr_read_type': {
        'check': ['clr', 'ccs', 'ul', 'reg']
    },
    'lr_input_format': {
        'check': ['fastq', 'pacbio_native']
    },
    'ss_library_fractions': {
        'check': [1, 2]
    },
    'nhr_assembler': {
        'check': ['flye', 'shasta', 'wtdbg', 'pereg', 'hifiasm'],
        'display': ['flye', 'shasta', 'wtdbg [wtdbg2]', 'pereg [Peregrine]', 'hifiasm']
    },
    'hap_assembler': {
        'check': ['flye', 'shasta', 'wtdbg', 'pereg', 'canu', 'hifiasm'],
        'display': ['flye', 'shasta', 'wtdbg [wtdbg2]', 'pereg [Peregrine]', 'canu', 'hifiasm']
    },
    'var_caller': {
        'check': ['freebayes', 'longshot', 'deepvar'],
        'display': ['freebayes', 'longshot', 'deepvar [DeepVariant]']
    }
}

option_help = col.OrderedDict([
    ('repo_folder', 'Full path to folder containing the Snakemake pipeline code '
                    '[Default: {}]'.format(option_defaults['repo_folder'])),
    ('exec_folder', 'Full path to base folder for pipeline execution '
                    '[Default: {}]'.format(option_defaults['exec_folder'])),
    ('sample_name', 'Name of the individual. Use only uppercase letters and numbers.'),
    ('super_population', 'Name or ID of super population. '
                         '[Default: {}]'.format(option_defaults['super_population'])),
    ('population', 'Name or ID of population [Default: {}]'.format(option_defaults['population'])),
    ('family_name', 'Name or ID of family [Default: {}]'.format(option_defaults['family_name'])),
    ('lr_folder', 'Full path to folder containing long-read input data. '
                  'All files in this folder are assumed to be long-read input files. '
                  'ATTENTION: remove any index or other auxiliary files from this folder. '
                  'Path must not be identical with the input data path for Strand-seq data.'),
    ('lr_project', 'Project name for long-read input data [Default: {}]'.format(option_defaults['lr_project'])),
    ('lr_seq_platform', 'Specify long-read sequencing '
                        'platform: {}'.format(option_constraints['lr_seq_platform']['display'])),
    ('lr_read_type', 'Specify type of long reads: {}'.format(option_constraints['lr_read_type']['check'])),
    ('lr_input_format', 'Specify input data format: {} (FASTQ must be '
                        'gzipped)'.format(option_constraints['lr_input_format']['check'])),
    ('ss_folder', 'Full path to folder containing Strand-seq input data (gzipped FASTQ). '
                  'All files in this folder are assumed to be Strand-seq read input files. '
                  'ATTENTION: remove any index or other auxiliary files from this folder. '
                  'Path must not be identical with the input path for long-read input data.'),
    ('ss_project', 'Project name for Strand-seq input data [Default: {}]'.format(option_defaults['ss_project'])),
    ('ss_seq_platform', 'Specify Strand-seq platform [Default: '
                        '{} (any Illumina)]'.format(option_defaults['ss_seq_platform'])),
    ('ss_read_type', 'Specify Strand-seq read type [Default: '
                     '{} (length N paired-end)]'.format(option_defaults['ss_read_type'])),
    ('ss_library_fractions', 'Specify number of different Strand-seq library fractions: '
                             '{}'.format(option_constraints['ss_library_fractions']['check'])),
    ('param_ver', 'Specify the pipeline parameter version (number) to be used in the analysis. '
                  '[Default: {}]'.format(option_defaults['param_ver'])),
    ('max_cpu', 'Number of CPU cores to use for tasks such as whole-genome assembly (should be max. available).'),
    ('high_cpu', 'Number of CPU cores to use for tasks such as whole-genome alignment of long reads.'),
    ('medium_cpu', 'Number of CPU cores to use for tasks such as chromosome-scale alignment of long reads.'),
    ('low_cpu', 'Number of CPU cores to use for tasks such as Strand-seq read alignment (should be 1 < X < 10).'),
    ('local_copy', 'Should the pipeline copy local input files (yes / no) '
                   '[Default: {}]?'.format(option_defaults['local_copy'])),
    ('singularity_module', 'What is the name of the ENV module that must be loaded '
                           'to run Singularity containers (needed for Peregrine and DeepVariant)? '
                           '[Default: <empty> ; do not use ENV module]'),
    ('nhr_assembler', 'Which assembler should be used to produce the initial non-haplotype resolved '
                      '("collapsed") assembly: {}'.format(option_constraints['nhr_assembler']['display'])),
    ('hap_assembler', 'Which assembler should be used to produce haploid '
                      'assemblies: {}'.format(option_constraints['hap_assembler']['display'])),
    ('var_caller', 'Which variant caller should be used: {}'.format(option_constraints['var_caller']['display']))
])

assert all(k in option_defaults.keys() for k in option_help.keys()), 'Key missing in option default/help'


load_configs = {
    'reference_sources': 'smk_config/ref_data/reference_data_sources.yml',
    'params': 'smk_config/params/smk_cfg_params_RV{parameter_version}.yml'
}

sample_targets_template = """
sample_targets_{sample_name}:
  - aliases:
      1: &long_reads {long_reads}
      2: &strandseq_reads {strandseq_reads}
  - defaults:
      hap_reads: *long_reads
      vc_reads: *long_reads
      sseq_reads: *strandseq_reads
      pol_reads: *long_reads
      hap_assm_mode: split
      hap:
        - h1-un
        - h2-un
  - target:
      nhr_assembler: {nhr_assembler}
      hap_assembler: {hap_assembler}
      var_caller: {var_caller}
      {pol_pass}
"""


def parse_command_line():

    parser = argparse.ArgumentParser(prog='autoconf.py', description=__desc__)

    parser.add_argument(
        '--accept-defaults',
        '-ad',
        action='store_true',
        default=False,
        dest='accept_defaults',
        help='Accept default values for all parameters if available.'
    )

    parser.add_argument(
        '--unique-id-length',
        '-uil',
        default=15,
        type=int,
        dest='unique_id_length',
        help='Use this many letters from the generated MD5 as unique ID for run and/or library.'
    )

    parser.add_argument(
        '--repository-folder',
        '-rf',
        default=option_defaults['repo_folder'],
        type=str,
        dest='repo_folder',
        help=option_help['repo_folder']
    )
    parser.add_argument(
        '--execution-folder',
        '-ef',
        default=option_defaults['exec_folder'],
        type=str,
        dest='exec_folder',
        help=option_help['exec_folder']
    )

    sample_group = parser.add_argument_group('Sample configuration')
    sample_group.add_argument(
        '--sample-name',
        '-sn',
        default=option_defaults['sample_name'],
        type=str,
        dest='sample_name',
        help=option_help['sample_name']
    )
    sample_group.add_argument(
        '--super-population',
        '-sp',
        default=option_defaults['super_population'],
        type=str,
        dest='super_population',
        help=option_help['super_population']
    )
    sample_group.add_argument(
        '--population',
        '-pop',
        default=option_defaults['population'],
        type=str,
        dest='population',
        help=option_help['population']
    )
    sample_group.add_argument(
        '--family-name',
        '-fam',
        default=option_defaults['family_name'],
        type=str,
        dest='family_name',
        help=option_help['family_name']
    )

    longread_group = parser.add_argument_group('Long-read input configuration')
    longread_group.add_argument(
        '--lr-data-folder',
        '-lrdf',
        default=option_defaults['lr_folder'],
        type=str,
        dest='lr_folder',
        help=option_help['lr_folder']
    )
    longread_group.add_argument(
        '--lr-project',
        '-lrp',
        default=option_defaults['lr_project'],
        type=str,
        dest='lr_project',
        help=option_help['lr_project']
    )
    longread_group.add_argument(
        '--lr-seq-platform',
        '-lrsp',
        default=option_defaults['lr_seq_platform'],
        type=str,
        choices=option_constraints['lr_seq_platform']['check'],
        dest='lr_seq_platform',
        help=option_help['lr_seq_platform']
    )
    longread_group.add_argument(
        '--lr-read-type',
        '-lrrt',
        default=option_defaults['lr_read_type'],
        type=str,
        choices=option_constraints['lr_read_type']['check'],
        dest='lr_read_type',
        help=option_help['lr_read_type']
    )
    longread_group.add_argument(
        '--lr-input-format',
        '-lrif',
        default=option_defaults['lr_input_format'],
        type=str,
        choices=option_constraints['lr_input_format']['check'],
        help=option_help['lr_input_format']
    )

    strandseq_group = parser.add_argument_group('Strand-seq input configuration')
    strandseq_group.add_argument(
        '--ss-data-folder',
        '-ssdf',
        default=option_defaults['ss_folder'],
        type=str,
        dest='ss_folder',
        help=option_help['ss_folder']
    )
    strandseq_group.add_argument(
        '--ss-project',
        '-ssp',
        default=option_defaults['ss_project'],
        type=str,
        dest='ss_project',
        help=option_help['ss_project']
    )
    strandseq_group.add_argument(
        '--ss-seq-platform',
        '-sssp',
        default=option_defaults['ss_seq_platform'],
        type=str,
        dest='ss_seq_platform',
        help=option_help['ss_seq_platform']
    )
    strandseq_group.add_argument(
        '--ss-read-type',
        '-ssrt',
        default=option_defaults['ss_read_type'],
        type=str,
        dest='ss_read_type',
        help=option_help['ss_read_type']
    )
    strandseq_group.add_argument(
        '--ss-library-fractions',
        '-sslf',
        default=option_defaults['ss_library_fractions'],
        type=int,
        choices=option_constraints['ss_library_fractions']['check'],
        dest='ss_library_fractions',
        help=option_help['ss_library_fractions']
    )

    param_group = parser.add_argument_group('Parameter configuration')
    param_group.add_argument(
        '--parameter-version',
        '-pv',
        default=option_defaults['param_ver'],
        type=int,
        dest='param_ver',
        help=option_help['param_ver']
    )

    runenv_group = parser.add_argument_group('Run environment configuration')
    runenv_group.add_argument(
        '--num-max-cpu',
        '-nmc',
        default=option_defaults['max_cpu'],
        type=int,
        dest='max_cpu',
        help=option_help['max_cpu']
    )
    runenv_group.add_argument(
        '--num-high-cpu',
        '-nhc',
        default=option_defaults['high_cpu'],
        type=int,
        dest='high_cpu',
        help=option_help['high_cpu']
    )
    runenv_group.add_argument(
        '--num-medium-cpu',
        '-ndc',
        default=option_defaults['medium_cpu'],
        type=int,
        dest='medium_cpu',
        help=option_help['medium_cpu']
    )
    runenv_group.add_argument(
        '--num-low-cpu',
        '-nlc',
        default=option_defaults['low_cpu'],
        type=int,
        dest='low_cpu',
        help=option_help['low_cpu']
    )
    runenv_group.add_argument(
        '--local-copy',
        '-lc',
        default=option_defaults['local_copy'],
        type=str,
        dest='local_copy',
        help=option_help['local_copy']
    )
    runenv_group.add_argument(
        '--singularity-module',
        '-sm',
        default=option_defaults['singularity_module'],
        type=str,
        dest='singularity_module',
        help=option_help['singularity_module']
    )

    build_target_group = parser.add_argument_group('Pipeline build targets')
    build_target_group.add_argument(
        '--nhr-assembler',
        '-na',
        default=option_defaults['nhr_assembler'],
        type=str,
        choices=option_constraints['nhr_assembler']['check'],
        dest='nhr_assembler',
        help=option_help['nhr_assembler']
    )
    build_target_group.add_argument(
        '--hap-assembler',
        '-ha',
        default=option_defaults['hap_assembler'],
        type=str,
        choices=option_constraints['hap_assembler']['check'],
        dest='hap_assembler',
        help=option_help['hap_assembler']
    )
    build_target_group.add_argument(
        '--var-caller',
        '-vc',
        default=option_defaults['var_caller'],
        type=str,
        choices=option_constraints['var_caller']['check'],
        dest='var_caller',
        help=option_help['var_caller']
    )

    args = parser.parse_args()
    return args


def collect_user_input(args):

    for option, help_text in option_help.items():
        option_default = option_defaults[option]
        option_set = getattr(args, option)
        if option_set is not None and args.accept_defaults:
            # Parameter has default value, may or may not have been changed
            # by the user
            # can skip to next
            continue
        elif option_set is not None and option_set != option_default:
            # do not accept default value...
            # user specified via command line
            # can skip to next
            continue
        else:
            # ask for user input
            pass
        print(help_text)
        value = input('--> ')
        if not value.strip() and option_default is None:
            raise ValueError('\n\nERROR: The parameter "{}" must not be empty\n\n'
                             'Parameter info: {}\n'.format(option, help_text))
        elif not value.strip():
            # user hit enter accepting default
            value = option_default
        else:
            pass
        setattr(args, option, value)

    return args


def sanitize_user_input(args):

    for option, help_text in option_help.items():
        value = getattr(args, option)

        if '_cpu' in option or option == 'ss_library_fractions':
            try:
                value = int(value)
            except ValueError:
                raise ValueError('\n\nParameter "{}" set to non-integer value "{}".\n\n'
                                 'Parameter info: {}\n\n'.format(option, value, help_text))

        constraints = option_constraints.get(option, {'check': []})['check']
        if constraints and value not in constraints:
            raise ValueError('\n\nParameter "{}" set to value "{}" - allowed values: {}\n\n'
                             'Parameter info: {}\n\n'.format(option, value, constraints, help_text))

        if option in ['repo_folder', 'lr_folder', 'ss_folder']:
            if not os.path.isdir(value):
                raise ValueError('\n\nThe folder "{}" (parameter "{}") does not exist.\n\n'
                                 'Parameter info: {}\n'.format(value, option, help_text))

        if option == 'local_copy':
            if not value or value.lower() in ['no', 'n', 'false', '0']:
                value = False
            else:
                value = True

        if option == 'singularity_module':
            if not bool(value):
                # empty string
                value = False

        setattr(args, option, value)

    if args.ss_folder.rstrip('/') == args.lr_folder.rstrip('/'):
        raise ValueError('\n\nERROR: The folder containing the Strand-seq input data '
                         'must not be identical with the long-read input '
                         'data folder: "{}" and "{}"\n\n'.format(args.ss_folder, args.lr_folder))
    return args


def collect_long_read_input(args):

    if args.lr_seq_platform.startswith('pb'):
        technology = 'pacbio'
    elif args.lr_seq_platform.startswith('ont'):
        technology = 'nanopore'
    else:
        raise ValueError('Unexpected long-read technology specified: {}'.format(args.lr_seq_platform))

    long_read_config = {
        'readset': '_'.join([args.sample_name, args.lr_project,
                             args.lr_seq_platform + '-' + args.lr_read_type]),
        'technology': technology,
        'data_type': args.lr_input_format,
        'load_type': None
    }

    input_files = sorted(filter(lambda x: os.path.isfile(os.path.join(args.lr_folder, x)), os.listdir(args.lr_folder)))

    if len(input_files) == 1:
        long_read_config['load_type'] = 'complete'
    elif len(input_files) > 1:
        long_read_config['load_type'] = 'parts'
    else:
        raise ValueError('No long-read input data found at path {}'.format(args.lr_folder))

    if long_read_config['data_type'] == 'fastq':
        file_ext = '.fastq.gz'
    elif long_read_config['data_type'] == 'pacbio_native':
        file_ext = '.pbn.bam'
    else:
        raise ValueError('Unexpected long-read data type: {}'.format(long_read_config['data_type']))

    if args.lr_read_type == 'clr' and args.lr_input_format == 'fastq':
        warnings.warn('Long-read data type was specified as FASTQ, but sequencing platform indicates '
                      'CLR reads. Polishing CLR assemblies requires PacBio-specific quality information '
                      'that is only contained in "pacbio-native" BAM files. You should absolutely be sure '
                      'of what you are doing.')

    link_files = []
    if long_read_config['load_type'] == 'complete':
        new_file_name = long_read_config['readset'] + '_1000' + file_ext
        new_file_path = os.path.join(args.exec_folder, 'autoconf_linked_data', new_file_name)
        old_file_path = os.path.join(args.lr_folder, input_files[0])
        link_files.append((old_file_path, new_file_path))
    else:
        new_file_top_path = os.path.join(args.exec_folder, 'autoconf_linked_data', long_read_config['readset'])
        for part_num, fname in enumerate(input_files, start=1):
            new_fname = long_read_config['readset'] + '.part{}'.format(part_num) + file_ext
            old_path = os.path.join(args.lr_folder, fname)
            link_files.append((old_path, os.path.join(new_file_top_path, new_fname)))

    return long_read_config, link_files


def find_mate_pairs(strandseq_files, unique_id_length):
    """
    Assume that second in pair is lexicographically
    closest to first in pair; just make a single confirmation
    (success or fail)
    :param strandseq_files: sorted list of Strand-seq input files
    :param unique_id_length:
    :return:
    """
    pairs = []
    run_ids = set()
    for pos in range(0, len(strandseq_files), 2):
        first_mate = strandseq_files[pos]
        second_mate = strandseq_files[pos+1]
        sm = difflib.SequenceMatcher(a=first_mate, b=second_mate)
        mobj = sm.find_longest_match(0, len(first_mate), 0, len(second_mate))
        match_size = mobj.size
        common_substring = first_mate[mobj.a : mobj.a+mobj.size]
        try:
            next_mate = strandseq_files[pos+2]
            sm.set_seq1(next_mate)
            mobj = sm.find_longest_match(0, len(next_mate), 0, len(second_mate))
            if match_size <= mobj.size:
                raise ValueError('Cannot pair mates based on lexicographical order: '
                                 '{} / {} vs {}'.format(first_mate, second_mate, next_mate))
        except IndexError:
            pass
        uniq_runid = hashlib.md5(common_substring.encode('utf-8')).hexdigest()[:unique_id_length]
        if uniq_runid in run_ids:
            raise ValueError('Created run ID is not unique: {} / {}'.format(uniq_runid, common_substring))
        run_ids.add(uniq_runid)
        pairs.append((first_mate, second_mate, common_substring, uniq_runid))
    return sorted(pairs, key=lambda x: x[2])


def match_library_fractions(mate_pairs, unique_id_length):
    """
    This might also handle files with unreasonable names
    :param mate_pairs:
    :param unique_id_length:
    :return:
    """
    last_fixed_id = mate_pairs[0][3]
    longest_match = 0
    longest_substring = ''
    matched_pair_ids = None

    matched_run_ids = set()
    generated_lib_ids = set()

    run_to_lib_id = dict()

    sm = difflib.SequenceMatcher(a=mate_pairs[0][2], b=mate_pairs[1][2])

    for pair_fix, pair_check in itertools.combinations(mate_pairs, 2):
        sub_fix, run_id_fix = pair_fix[2], pair_fix[3]
        sub_check, run_id_check = pair_check[2], pair_check[3]
        if run_id_fix != last_fixed_id:
            if matched_pair_ids is not None:
                uniq_lib_id = hashlib.md5(longest_substring.encode('utf-8')).hexdigest()[:unique_id_length]
                if uniq_lib_id in generated_lib_ids:
                    raise ValueError('Iter -> created lib ID is not unique: {} / {}'.format(uniq_lib_id, longest_substring))
                # associate run id to lib id
                run_to_lib_id[matched_pair_ids[0]] = uniq_lib_id
                run_to_lib_id[matched_pair_ids[1]] = uniq_lib_id

                matched_run_ids.add(matched_pair_ids[0])
                matched_run_ids.add(matched_pair_ids[1])

                generated_lib_ids.add(uniq_lib_id)

            # reset everything
            matched_pair_ids = None
            longest_match = 0
            longest_substring = ''

            last_fixed_id = run_id_fix
            sm = difflib.SequenceMatcher(a=sub_check, b=sub_fix)
        else:
            sm.set_seq1(sub_check)
        if run_id_fix in matched_run_ids or run_id_check in matched_run_ids:
            continue
        mobj = sm.find_longest_match(0, len(sub_check), 0, len(sub_fix))
        if longest_match < mobj.size:
            matched_pair_ids = run_id_fix, run_id_check
            longest_match = mobj.size
            longest_substring = sub_fix[mobj.b:mobj.b + mobj.size]

    uniq_lib_id = hashlib.md5(longest_substring.encode('utf-8')).hexdigest()[:unique_id_length]
    if uniq_lib_id in generated_lib_ids:
        raise ValueError('Last iter -> created lib ID is not unique: {} / {}'.format(uniq_lib_id, longest_substring))
    # associate run id to lib id
    run_to_lib_id[matched_pair_ids[0]] = uniq_lib_id
    run_to_lib_id[matched_pair_ids[1]] = uniq_lib_id

    matched_run_ids.add(matched_pair_ids[0])
    matched_run_ids.add(matched_pair_ids[1])

    assert len(mate_pairs) == len(matched_run_ids), \
        'Missing matched library fractions: expect {} / have {}'.format(len(mate_pairs), len(matched_run_ids))

    return run_to_lib_id


def compile_new_strandseq_path(base_folder, readset, lib_id, run_id, mate_num):

    if lib_id is not None:
        matched_lib = 'L' + str(lib_id[run_id])
        fname = '_'.join([readset, matched_lib, 'R' + str(run_id), mate_num]) + '.fastq.gz'
    else:
        fname = '_'.join([readset, 'R' + str(run_id), mate_num]) + '.fastq.gz'
    full_path = os.path.join(base_folder, 'autoconf_linked_data', readset + '_sseq', fname)
    return full_path


def collect_strandseq_input(args):

    lib_frac = {1: 'one', 2: 'two'}[args.ss_library_fractions]
    strandseq_config = {
        'readset': '_'.join([args.sample_name, args.ss_project,
                             args.ss_seq_platform + '-' + args.ss_read_type, 'sseq']),
        'library_fractions': lib_frac
    }

    strandseq_files = sorted(filter(lambda x: os.path.isfile(os.path.join(args.ss_folder, x)), os.listdir(args.ss_folder)))

    pairs = find_mate_pairs(strandseq_files, args.unique_id_length)
    lib_ids = None
    if args.ss_library_fractions > 1:
        lib_ids = match_library_fractions(pairs, args.unique_id_length)

    link_files = []

    file_readset = strandseq_config['readset'].rsplit('_', 1)[0]
    for mate1, mate2, _, run_id in pairs:
        old_path1 = os.path.join(args.ss_folder, mate1)
        new_path1 = compile_new_strandseq_path(args.exec_folder, file_readset, lib_ids, run_id, '1')
        link_files.append((old_path1, new_path1))

        old_path2 = os.path.join(args.ss_folder, mate2)
        new_path2 = compile_new_strandseq_path(args.exec_folder, file_readset, lib_ids, run_id, '2')
        link_files.append((old_path2, new_path2))

    return strandseq_config, link_files


def define_default_polishing(args):

    if args.lr_read_type == 'clr':
        pol_pass = 'pol_pass: arrow-p1'
    elif args.lr_read_type == 'ccs':
        pol_pass = 'pol_pass: racon-p2'
    else:
        pol_pass = ''
    return pol_pass


def generate_sample_description(args, lr_readset, sseq_readset):

    sample_desc = {
        'sample_description_{}'.format(args.sample_name): dict([
            ('individual', args.sample_name),
            ('super_population', args.super_population),
            ('population', args.population),
            ('family', args.family_name),
            ('data_sources', [
                {'long_reads': lr_readset},
                {'strandseq': sseq_readset}
            ])
        ])}
    return sample_desc


def generate_sample_targets(args, lr_readset, sseq_readset):

    lr_readset_name = lr_readset + '_1000'
    sseq_readset_name = sseq_readset
    polishing = define_default_polishing(args)

    sample_targets = sample_targets_template.format(**{
        'sample_name': args.sample_name,
        'long_reads': lr_readset_name,
        'strandseq_reads': sseq_readset_name,
        'nhr_assembler': args.nhr_assembler,
        'hap_assembler': args.hap_assembler,
        'var_caller': args.var_caller,
        'pol_pass': polishing
    })

    return sample_targets.strip()


def generate_run_env(args):

    run_env = {
        'env_module_singularity': args.singularity_module,
        'force_local_copy': args.local_copy,
        'num_cpu_max': args.max_cpu,
        'num_cpu_high': args.high_cpu,
        'num_cpu_medium': args.medium_cpu,
        'num_cpu_low': args.low_cpu
    }
    return run_env


def generate_data_source(args):

    bam_extension = 'pbn.bam' if args.lr_input_format == 'pacbio_native' else 'bam'

    data_source = {
        'data_source_{}'.format(args.sample_name): {
            'output': '{}_local_data.json'.format(args.sample_name),
            'server': 'localhost',
            'data_source': os.path.join(args.exec_folder, 'autoconf_linked_data'),
            'collect_files': [
                'fastq.gz',
                bam_extension
            ],
            'sort_into': [
                'fastq',
                'bam'
            ],
            'assume_correct_filenames': True
        }
    }
    return data_source


def main():
    args = parse_command_line()
    args = collect_user_input(args)
    args = sanitize_user_input(args)

    lr_readset_config, lr_link_files = collect_long_read_input(args)
    sseq_readset_config, sts_link_files = collect_strandseq_input(args)

    config_path = os.path.join(args.exec_folder, 'autoconf_config', 'run_assembly.yml')
    os.makedirs(os.path.dirname(config_path), exist_ok=True)

    with open(config_path, 'w') as config_dump:

        sample_desc = generate_sample_description(args, lr_readset_config, sseq_readset_config)
        _ = config_dump.write('# SECTION: SAMPLE & READSET DESCRIPTION\n')
        _ = config_dump.write(yaml.dump(sample_desc, default_flow_style=False, sort_keys=False))
        _ = config_dump.write('\n\n')

        sample_targets = generate_sample_targets(args, lr_readset_config['readset'], sseq_readset_config['readset'])
        _ = config_dump.write('# SECTION: SAMPLE TARGETS / SNAKEMAKE WILDCARDS\n')
        _ = config_dump.write(sample_targets)
        _ = config_dump.write('\n\n')

        run_env = generate_run_env(args)
        _ = config_dump.write('# SECTION: RUN ENVIRONMENT\n')
        _ = config_dump.write(yaml.dump(run_env, default_flow_style=False, sort_keys=False))
        _ = config_dump.write('\n\n')

        data_source = generate_data_source(args)
        _ = config_dump.write('# SECTION: DATA SOURCE\n')
        _ = config_dump.write(yaml.dump(data_source, default_flow_style=False, sort_keys=False))
        _ = config_dump.write('\n\n')

        _ = config_dump.write('# SECTION: PIPELINE PARAMETERS\n')
        with open(os.path.join(args.repo_folder, load_configs['params'].format(parameter_version=args.param_ver)), 'r') as cfg:
            _ = config_dump.write(cfg.read())
            _ = config_dump.write('\n\n')

        _ = config_dump.write('# SECTION: REFERENCE DATA\n')
        with open(os.path.join(args.repo_folder, load_configs['reference_sources']), 'r') as cfg:
            _ = config_dump.write(cfg.read())
            _ = config_dump.write('\n')

    # perform load check
    with open(config_path, 'r') as cfg:
        _ = yaml.load(cfg, Loader=yaml.SafeLoader)

    for old_path, new_path in itertools.chain(lr_link_files, sts_link_files):
        os.makedirs(os.path.dirname(new_path), exist_ok=True)
        os.symlink(old_path, new_path)

    return 0


if __name__ == '__main__':
    main()
