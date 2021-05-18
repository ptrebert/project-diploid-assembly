
import os
import sys
import collections as col

localrules: master_scrape_data_sources


def read_configured_data_sources():
    """
    :return:
    """

    COMMAND_SCAN_PATH_CALL = '{{script_exec}} --debug --output {{output}} '
    COMMAND_SCAN_PATH_CALL += '--server {server} --data-source {data_source} '
    COMMAND_SCAN_PATH_CALL += '--collect-files {collect_files} --sort-into {sort_into} '
    COMMAND_SCAN_PATH_CALL += '{file_infix} {file_suffix} '
    COMMAND_SCAN_PATH_CALL += '{local_path_suffix} {fix_tech} '
    COMMAND_SCAN_PATH_CALL += '{assume_pacbio_native} {assume_clr_subreads} '
    COMMAND_SCAN_PATH_CALL += '{assume_paired_reads} {assume_correct_filenames} '
    COMMAND_SCAN_PATH_CALL += '&> {{log}}'

    default_parameters = {
        'file_infix': '',
        'file_suffix': '',
        'local_path_suffix': '',
        'fix_tech': '',
        'assume_pacbio_native': '',
        'assume_clr_subreads': '',
        'assume_paired_reads': '',
        'assume_correct_filenames': ''
    }

    script_path = find_script_path('scan_remote_path.py')

    source_names = []
    source_outputs = []
    source_syscalls = []

    for key, values in config.items():
        if not key.startswith('data_source'):
            continue
        source_name = key.split('_', 2)[-1]

        data_source_output = None
        data_source_config = dict(default_parameters)

        for parameter, setting in values.items():
            if parameter == 'output':
                if not setting.startswith('input/data_sources'):
                    data_source_output = os.path.join('input/data_sources', setting)
                else:
                    data_source_output = setting
            elif parameter.startswith('assume'):
                bool_switch = ''.join(['--', parameter.replace('_', '-')])
                if setting:
                    data_source_config[parameter] = bool_switch
                else:
                    data_source_config[parameter] = ''
            elif any([parameter.startswith(x) for x in ['name', 'comment']]):
                source_name += ' ({})'.format(setting)
            else:
                if isinstance(setting, list):
                    setting = ' '.join(setting)
                if parameter in default_parameters:
                    data_source_config[parameter] = ''.join(['--', parameter.replace('_', '-')]) + ' ' + setting
                else:
                    data_source_config[parameter] = setting

        if data_source_output is None:
            sys.stderr.write('\nERROR: no output ("output: filename") file defined for data source: {}'.format(key))
            raise ValueError('No output file defined for data source: {}'.format(key))

        try:
            source_call = COMMAND_SCAN_PATH_CALL.format(**data_source_config)
        except KeyError as kerr:
            sys.stderr.write('\nERROR: cannot parse data source specification: {} / {}\n'.format(key, values))
            raise kerr
        else:
            source_names.append(source_name)
            source_outputs.append(data_source_output)
            log_file = os.path.join('log', data_source_output.replace('.json', '.log'))
            source_call = source_call.format(**{'script_exec': script_path,
                                                'output': data_source_output,
                                                'log': log_file})
            source_syscalls.append(source_call)

    if len(source_outputs) != len(set(source_outputs)):
        sys.stderr.write('\nERROR: at least two data sources specify the same output file: {}\n'.format(sorted(source_outputs)))
        raise ValueError('Duplicate data source output file.')

    if not source_names:
        if bool(config.get('show_warnings', False)):
            sys.stderr.write('\nWARNING: no data sources configured\n')

    if not len(source_names) == len(source_outputs) == len(source_syscalls):
        raise RuntimeError('Length mismatch for data sources: \n'
                           ' => {} \n => {} \n => {}'.format(source_names, source_outputs, source_syscalls))

    return source_names, source_outputs, source_syscalls


def collect_input_data_files(top_path, sample, file_ext):
    """
    """
    selected_files = []
    for root, dirs, files in os.walk(top_path, followlinks=False):
        if 'ignore' in root:
            continue
        level_files = [f for f in files if sample in f and f.endswith(file_ext)]
        if not level_files:
            level_files = [f for f in files if f.endswith(file_ext)]
        level_files = [(f, os.path.join(root, f)) for f in level_files]
        selected_files.extend(level_files)
    selected_files = sorted(set(selected_files))
    if len(selected_files) < 1:
        raise ValueError('No data files found: {} / {} / {}'.format(top_path, sample, file_ext))
    return selected_files


def check_readset_name(readset_type, readset_name, sample_name):
    """
    """
    if not readset_name.startswith(sample_name):
        raise ValueError('Readset name has to start with sample name: {} / {}'.format(readset_name, sample_name))
    if readset_type == 'long_reads':
        try:
            sample, project, platform = readset_name.split('_')
        except ValueError:
            raise ValueError('Readset name does not follow naming convention: SAMPLE_PROJECT_SEQPLATFORM')
    
    elif readset_type in ['strandseq', 'short_reads']:
        try:
            sample, project, platform, suffix = readset_name.split('_')
        except ValueError:
            expected_suffix = 'sseq' if readset_type == 'strandseq' else 'short'
            raise ValueError('Readset name does not follow naming convention: SAMPLE_PROJECT_SEQPLATFORM_{}'.format(expected_suffix))
    else:
        raise ValueError('Unknown type of readset: {} / {}'.format(readset_type, readset_name))


def determine_data_source_file_extension(readset_spec):
    """
    """

    has_data_type = 'data_type' in readset_spec

    if 'data_source_filter' in readset_spec:
        file_ext = readset_spec['data_source_filter'].strip().strip('"')
        if has_data_type and readset_spec['data_type'] == 'pacbio_native':
            sort_folder = 'input/bam'
            pgas_ext = '.pbn.bam'
        elif (not has_data_type) or readset_spec['data_type'] == 'fastq':
            sort_folder = 'input/fastq'
            pgas_ext = '.fastq.gz'
        else:
            raise ValueError('Cannot process readset data type: {} / {}'.format(readset_spec['readset'], readset_spec['data_type']))
    elif has_data_type and readset_spec['data_type'] == 'pacbio_native':
        file_ext = '.bam'
        sort_folder = 'input/bam'
        pgas_ext = '.pbn.bam'
    elif (not has_data_type) or readset_spec['data_type'] == 'fastq':
        file_ext = '.fastq.gz'
        sort_folder = 'input/fastq'
        pgas_ext = '.fastq.gz'
    else:
        raise ValueError('Cannot determine file extension to filter data source: {}'.format(readset_spec))
    return file_ext, sort_folder, pgas_ext


def handle_long_read_data_source(readset_spec, sample_name):
    """
    """
    file_ext, sort_folder, pgas_ext = determine_data_source_file_extension(readset_spec)
    input_data_files = collect_input_data_files(
        readset_spec['data_source_folder'],
        sample_name,
        file_ext
    )
    readset_name = readset_spec['readset']
    data_sources = col.OrderedDict()
    if len(input_data_files) == 1:
        if not readset_spec['load_type'] == 'complete':
            raise ValueError('Only one input file for readset {}, but specified as coming in "parts" '
                            '- incomplete input data?'.format(readset_name))
        req_key = os.path.join(sort_folder, readset_name + '_1000')
        data_sources[req_key] = {
                'local_path': os.path.join(sort_folder, readset_name + '_1000' + pgas_ext),
                'remote_path': input_data_files[0][1]
            }
    else:
        for pos, (fname, remote_path) in enumerate(input_data_files, start=1):
            part_file = readset_name + '.part' + str(pos)
            req_key = os.path.join(sort_folder, readset_name, part_file)
            data_sources[req_key] = {
                'local_path': os.path.join(sort_folder, readset_name, part_file + pgas_ext),
                'remote_path': remote_path
            }
    return data_sources


def compile_data_source_entry(sort_folder, readset_name, readset_prefix, fraction, library_id, run_id):
    """
    """
    if fraction is None:
        read1_file_prefix = '_'.join([readset_prefix, library_id, '1'])
        read1_key = os.path.join(sort_folder, readset_name, read1_file_prefix)

        read2_file_prefix = '_'.join([readset_prefix, library_id, '2'])
        read2_key = os.path.join(sort_folder, readset_name, read2_file_prefix)

    else:
        read1_file_prefix = '_'.join([readset_prefix + '-' + fraction, library_id, run_id, '1'])
        read1_key = os.path.join(sort_folder, readset_name, read1_file_prefix)

        read2_file_prefix = '_'.join([readset_prefix + '-' + fraction, library_id, run_id, '2'])
        read2_key = os.path.join(sort_folder, readset_name, read2_file_prefix)

    return read1_key, read2_key


def group_by_two_iterator(file_list, readset_name, sort_folder, pgas_ext):
    """
    """
    import hashlib

    if len(file_list) % 2 != 0:
        raise ValueError('Number of data files cannot be grouped by {} ({})'.format(2, len(file_list)))

    readset_prefix = readset_name.rsplit('_', 1)[0]

    for i in range(0, len(file_list), 2):
        (read1_fname, read1_fpath), (read2_fname, read2_fpath) = file_list[i:i+2]
        library_id = ''.join([read1_fname, read2_fname])
        library_id = hashlib.md5(library_id.encode('utf-8')).hexdigest()

        read1_key, read2_key = compile_data_source_entry(
            sort_folder,
            readset_name,
            readset_prefix,
            None,
            library_id,
            None,
        )

        read_info = {
            read1_key: {
                'local_path': read1_key + pgas_ext,
                'remote_path': read1_fpath
            },
            read2_key: {
                'local_path': read2_key + pgas_ext,
                'remote_path': read2_fpath
            }
        }
        
        yield read_info

    return


def group_by_four_iterator(file_list, readset_name, sort_folder, pgas_ext):
    """
    """
    import hashlib

    if len(file_list) % 4 != 0:
        raise ValueError('Number of data files cannot be grouped by {} ({})'.format(4, len(file_list)))

    # - This is different from "grouping by two" because we need to account
    # here for the two different sequencing fraction per library
    # - This cannot fail because it has already been checked before
    sample, project, platform, suffix = readset_name.split('_')
    if '-' not in platform:
        raise ValueError('For Strand-seq samples with two sequencing fractions per library, '
                        'the platform string ({}) needs to be separable by "minus" to encode '
                        'the two fraction (i.e. SOMETHING-frac1 and SOMETHING-frac2'.format(platform))
    sequencer = platform.split('-')[0]
    readset_prefix = '_'.join([sample, project, sequencer])

    for i in range(0, len(file_list), 4):
        read_file_names = [t[0] for t in file_list[i:i+4]]
        read_file_paths = [t[1] for t in file_list[i:i+4]]
        
        library_id = ''.join(read_file_names)
        library_id = hashlib.md5(library_id.encode('utf-8')).hexdigest()

        frac1_run_id = ''.join(read_file_names[:2])
        frac1_run_id = hashlib.md5(frac1_run_id.encode('utf-8')).hexdigest()

        frac1_read1_key, frac1_read2_key = compile_data_source_entry(
            sort_folder,
            readset_name,
            readset_prefix,
            'frac1',
            library_id,
            frac1_run_id,
        )

        frac2_run_id = ''.join(read_file_names[2:])
        frac2_run_id = hashlib.md5(frac2_run_id.encode('utf-8')).hexdigest()

        frac2_read1_key, frac2_read2_key = compile_data_source_entry(
            sort_folder,
            readset_name,
            readset_prefix,
            'frac2',
            library_id,
            frac2_run_id,
        )

        read_info = {
            frac1_read1_key: {
                'local_path': frac1_read1_key + pgas_ext,
                'remote_path': read_file_paths[0]
            },
            frac1_read2_key: {
                'local_path': frac1_read2_key + pgas_ext,
                'remote_path': read_file_paths[1]
            },
            frac2_read1_key: {
                'local_path': frac2_read1_key + pgas_ext,
                'remote_path': read_file_paths[2]
            },
            frac2_read2_key: {
                'local_path': frac2_read2_key + pgas_ext,
                'remote_path': read_file_paths[3]
            }
        }
              
        yield read_info

    return


def handle_strandseq_data_source(readset_spec, sample_name):
    """
    """
    file_ext, sort_folder, pgas_ext = determine_data_source_file_extension(readset_spec)
    input_data_files = collect_input_data_files(
        readset_spec['data_source_folder'],
        sample_name,
        file_ext
    )
    readset_name = readset_spec['readset']
    data_sources = col.OrderedDict()
    num_lib_fractions = readset_spec['library_fractions']
    # TODO: with a little more abstraction, the below iterators
    # could be realized as a single function
    if num_lib_fractions == 'one':
        for read_info in group_by_two_iterator(input_data_files, readset_name, sort_folder, pgas_ext):
            data_sources.update(read_info)

    elif num_lib_fractions == 'two':
        for read_info in group_by_four_iterator(input_data_files, readset_name, sort_folder, pgas_ext):
            data_sources.update(read_info)
    else:
        raise ValueError('Cannot process number of Strand-seq library fractions: {}'.format(readset_spec))
    return data_sources


def handle_short_read_data_source(readset_spec, sample_name):
    raise NotImplementedError('TODO: implement short read input data')


def read_local_data_sources():
    """
    Create a quick and dirty JSON annotation for locally available data sources.
    This relies on a proper user annotation of readsets in the sample config.
    """

    readset_handlers = {
        'long_reads': handle_long_read_data_source,
        'strandseq': handle_strandseq_data_source,
        'short_reads': handle_short_read_data_source,
    }

    show_warnings = bool(config.get('show_warnings', False))

    formatted_data_sources = dict()

    for key, values in config.items():
        if not key.startswith('sample_description'):
            continue
        sample = key.split('_')[-1]

        for ds in values['data_sources']:
            for readset_type, readset_spec in ds.items():
                readset_name = readset_spec['readset']
                if 'data_source_folder' not in readset_spec:
                    if show_warnings:
                        sys.stderr.write('\nWarning: no data source folder for readset {} - skipping\n'.format(readset_name))
                    continue
                if not os.path.isdir(readset_spec['data_source_folder']):
                    if show_warnings:
                        sys.stderr.write('\nWarning: skipping over non-existing data source folder: {}\n'.format(readset_spec['data_source_folder']))
                    continue
                if readset_name in formatted_data_sources:
                    raise ValueError('Duplicate data source entry for readset {}'.format(readset_name))
                check_readset_name(readset_type, readset_name, sample)
                readset_sources = readset_handlers[readset_type](readset_spec, sample)
                formatted_data_sources[readset_name] = readset_sources
    
    return formatted_data_sources


USE_LEGACY_DATA_SCRAPING = bool(config.get('use_legacy_data_scraping', False))

if USE_LEGACY_DATA_SCRAPING:

    DATA_SOURCE_NAMES, DATA_SOURCE_OUTPUTS, DATA_SOURCE_CALLS = read_configured_data_sources()

    DATA_SOURCE_TO_CALL = dict((os.path.basename(source).rsplit('.', 1)[0], call) for source, call in zip(DATA_SOURCE_OUTPUTS, DATA_SOURCE_CALLS))

    FORMATTED_DATA_SOURCES = None

else:
    FORMATTED_DATA_SOURCES = read_local_data_sources()
    DATA_SOURCE_FILES = sorted(list(FORMATTED_DATA_SOURCES.keys()))
    DATA_SOURCE_OUTPUTS = ['input/data_sources/{}.json'.format(fn) for fn in DATA_SOURCE_FILES]
    DATA_SOURCE_TO_CALL = dict()


def get_strandseq_library_info(sseq_reads):
    """
    Helper function to extract needed information about Strand-seq
    libraries - used in collect strand-seq merge files and
    collect strand-seq alignments
    """
    import os
    import sys

    debug = bool(config.get('show_debug_messages', False))

    if debug:
        sys.stderr.write('Collecting Strand-seq library info for: {}\n'.format(sseq_reads))

    if FORMATTED_DATA_SOURCES is None:
        raise RuntimeError('No data sources formatted - is legacy data scraping set?')
    try:
        sseq_source_info = FORMATTED_DATA_SOURCES[sseq_reads]
    except KeyError:
        raise KeyError('No data sources for SSEQ readset found: {}'.format(sseq_reads))

    if sseq_reads in CONSTRAINT_STRANDSEQ_MONOFRACTION_SAMPLES:
        lib_id_pos = -2
        if debug:
            sys.stderr.write('Strand-seq readset is a mono-fraction sample\n'.format(sseq_reads))
    else:
        lib_id_pos = -3
        if debug:
            sys.stderr.write('Strand-seq readset is a dual-fraction sample\n'.format(sseq_reads))

    sseq_lib_infos = set()
    for sseq_lib, sseq_lib_source in sseq_source_info.items():
        # e.g.: HG00733_1kg_il25k-frac1_ab50c4c67063fdb907a4f49512a34f0e_c52ee9c925d8c7d8d7e534411e21495b_2
        lib_pair_name = os.path.basename(sseq_lib)
        # e.g.: HG00733_1kg_il25k-frac1_ab50c4c67063fdb907a4f49512a34f0e_c52ee9c925d8c7d8d7e534411e21495b
        lib_name = lib_pair_name.rsplit('_', 1)[0]
        # e.g.: ab50c4c67063fdb907a4f49512a34f0e
        lib_id = lib_pair_name.split('_')[lib_id_pos]
        # e.g.: ERR1295562_2.fastq.gz
        remote_source = os.path.basename(sseq_lib_source['remote_path'])
        
        sseq_lib_infos.add(
            (
                lib_name,
                lib_id,  # only needed for merging step for dual-fraction SSEQ samples
                remote_source  # this will be removed before returning
            )
        )
    
    if debug:
        sys.stderr.write('Collected info for {} libraries ({} files)\n'.format(len(sseq_lib_infos) // 2, len(sseq_lib_infos)))

    # if automatic library QC is set to yes for this sample,
    # need to subtract all excluded libraries if QC is complete
    # to avoid pipeline failures after restart
    exclude_libs = []
    if sseq_reads in CONSTRAINT_STRANDSEQ_LIBQC_SAMPLES:
        # check if exclude list already exists
        exclude_file = 'output/sseq_qc/{}.exclude.txt'.format(sseq_reads)
        if os.path.isfile(exclude_file):
            with open(exclude_file, 'r') as listing:
                exclude_libs.extend(listing.read().strip().split())
            if debug:
                sys.stderr.write('Loaded {} SSEQ library exclude hints from file: {}\n'.format(len(exclude_libs), exclude_file))
        else:
            if debug:
                sys.stderr.write('SSEQ exclude file does not exist yet - may cause job failures: {}\n'.format(exclude_file))

    sseq_libs = set()
    sseq_lib_ids = set()
    for lib_name, lib_id, remote_source in sorted(sseq_lib_infos):
        if any(x in remote_source for x in exclude_libs):
            if debug:
                sys.stderr.write('Excluding library {} for SSEQ readset {}\n'.format(lib_name, sseq_reads))
            continue
        sseq_libs.add(lib_name)
        sseq_lib_ids.add(lib_id)        
    
    if debug:
        sys.stderr.write('Returning {} selected libraries for readset {}\n'.format(len(sseq_libs), sseq_reads))

    if sseq_reads in CONSTRAINT_STRANDSEQ_MONOFRACTION_SAMPLES:
        assert len(sseq_libs) == len(sseq_lib_ids), \
            'Mismatch library names ({}) and IDs ({}) - name extract failed?'.format(len(sseq_libs), len(sseq_lib_ids))
    else:
        # for dual-fraction samples, there are 2 pairs of read files
        # that together make one library in the context of PGAS
        assert len(sseq_libs) == (len(sseq_lib_ids) * 2), \
            'Mismatch library names ({}) and IDs ({}) - name extract failed?'.format(len(sseq_libs), len(sseq_lib_ids))

    return sseq_libs, sseq_lib_ids


def count_number_of_input_parts(readset):
    """
    Helper function to get number of parts for
    (usually) long read input reads
    """
    if FORMATTED_DATA_SOURCES is None:
        raise RuntimeError('No data sources formatted - is legacy data scraping set?')
    try:
        read_source_info = FORMATTED_DATA_SOURCES[readset]
    except KeyError:
        raise KeyError('No data sources for readset found: {}'.format(readset))

    assert readset in CONSTRAINT_PARTS_FASTQ_INPUT_SAMPLES, 'This is not a PARTS/FASTQ input sample: {}'.format(readset)

    return len(read_source_info.values())



rule master_scrape_data_sources:
    input:
        DATA_SOURCE_OUTPUTS


if USE_LEGACY_DATA_SCRAPING:
    rule scrape_data_source:
        """
        2020-02-06
        This was originally an anon rule iterating through the above lists.
        This worked locally, but resulted in wrong system calls matched to certain
        output file names in a cluster environment (maybe an issue with the job
        submission script for anon rules?). Don't wait for snakemake fix, work around...
        """
        output:
            'input/data_sources/{data_source}.json'
        log:
            'log/input/data_sources/{data_source}.log'
        message: 'Processing data source: {output}'
        conda:
            '../environment/conda/conda_pyscript.yml'
        params:
            scrape_call = lambda wildcards: DATA_SOURCE_TO_CALL[wildcards.data_source]
        shell:
            '{params.scrape_call}'

else:
    rule scan_local_data_sources:
        """
        Check all sample configs for data paths and build a quick and dirty
        data source json from the gathered info.
        """
        output:
            'input/data_sources/{data_source}.json'
        log:
            'log/input/data_sources/{data_source}.log'
        message: 'Processing local data sources for PGAS run: {wildcards.data_source}'
        run:
            import json

            local_data_sources = FORMATTED_DATA_SOURCES

            if not local_data_sources:
                raise ValueError('\nERROR: no local data sources could be configured. Check your sample configuration.\n')

            for readset, input_sources in local_data_sources.items():
                if readset != wildcards.data_source:
                    continue
                with open(output[0], 'w') as dump:
                    json.dump(
                        input_sources,
                        dump,
                        ensure_ascii=True,
                        indent=1,
                        sort_keys=False
                    )            
