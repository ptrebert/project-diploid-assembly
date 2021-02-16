
localrules: master_handle_data_download


rule master_handle_data_download:
    input:
        rules.master_link_data_sources.output,


def create_request_files_from_json(json_dumps, request_path, blacklist, logfile):
    """
    Process output of scan_remote_path scraping script
    :param json_dumps:
    :param request_path:
    :param logfile:
    :return:
    """
    import os
    import errno
    import json as json
    import collections as col

    complete_sources = dict()
    check_key_dups = col.Counter()

    for annotation in json_dumps:
        _ = logfile.write('Reading json dump {}\n'.format(annotation))
        with open(annotation, 'r') as paths:
            obj = json.load(paths)
            check_key_dups.update(list(obj.keys()))
            complete_sources.update(obj)

    # add all manually annotated data sources
    for k, v in config.items():
        if not k.startswith('sample_data_sources'):
            continue
        manual_annotation = config[k]
        check_key_dups.update(list(manual_annotation.keys()))
        complete_sources.update(manual_annotation)

    count_keys = check_key_dups.most_common()
    if len(count_keys) > 1:
        if count_keys[0][1] > 1:
            dup_keys = [t for t in count_keys if t[1] > 1]
            raise ValueError('Duplicate keys in data sources annotation: {}'.format(dup_keys))

    path_sample_prefix = os.path.split(request_path)[0]
    requests_created = False

    for key, file_infos in complete_sources.items():
        if not key.startswith(path_sample_prefix):
            continue
        if not requests_created:
            _ = logfile.write('Creating request path {}\n'.format(request_path))
            os.makedirs(request_path, exist_ok=True)
            requests_created = True

        _ = logfile.write('Processing file key {}\n'.format(key))
        file_prefix = os.path.split(key)[1]
        req_file_path = os.path.join(request_path, file_prefix + '.request')
        try:
            md5 = file_infos['md5']
        except KeyError:
            md5 = 'no_md5'

        if blacklist and any([hint in file_infos['remote_path'] for hint in blacklist]):
            _ = logfile.write('File blacklisted: {}\n'.format(file_infos['local_path']))
            try:
                os.remove(req_file_path)
            except (OSError, IOError) as err:
                if err.errno != errno.ENOENT:
                    raise RuntimeError('Could not remove blacklisted request file {}: {}'.format(req_file_path, str(err)))
            try:
                os.remove(file_infos['local_path'])
            except (OSError, IOError) as err:
                if err.errno != errno.ENOENT:
                    raise RuntimeError('Could not remove blacklisted data file {}: {}'.format(file_infos['local_path'], str(err)))
            continue

        _ = logfile.write('Creating request file at {}\n'.format(req_file_path))
        with open(req_file_path, 'w') as req_file:
            _ = req_file.write(file_infos['remote_path'] + '\n')
            _ = req_file.write(file_infos['local_path'] + '\n')
            _ = req_file.write(md5 + '\n')
    return requests_created


def create_request_files_from_tsv(tsv_files, request_path, blacklist, logfile):
    """
    Process output of downloader / ena-file-report
    Note that this relies on the fact that the relevant
    tsv file is named after the readset (and not the bioproject),
    thus giving an exact match for the request path
    :param tsv_files:
    :param request_path:
    :param logfile:
    :return:
    """
    import os
    import errno
    import csv

    requests_created = False

    long_read_parts = 0

    mate1_indicators = [(1, x) for x in ['.R1.fastq.gz', '_R1.fastq.gz', '_1.fastq.gz']]
    mate2_indicators = [(2, x) for x in ['.R2.fastq.gz', '_R2.fastq.gz', '_2.fastq.gz']]

    for tsv in tsv_files:
        assert tsv.endswith('.metadata.tsv'), 'Unexpected file name: {}'.format(tsv)
        filename = os.path.basename(tsv).split('.')[0]  # strip .metadata.tsv
        if filename not in request_path:
            _ = logfile.write('Skipping over TSV file {}\n'.format(tsv))
            continue
        individual = filename.split('_')[0]
        sample_annotator = None
        with open(tsv, 'r', newline='') as table:
            rows = csv.DictReader(table, delimiter='\t')
            for row in rows:
                if individual not in row['sample_alias']:
                    continue
                if sample_annotator is None:
                    bioproject = row.get('study_accession', False)
                    if not bioproject:
                        raise ValueError('Cannot identify bioproject/study accession in row: {}'.format(row))
                    sample_annotator = get_bioproject_sample_annotator(bioproject)
                label, tech_type = sample_annotator(row, individual)
                if label is None:
                    continue
                if tech_type == 'long':
                    long_read_parts += 1
                _ = logfile.write('Processing data for sample {}\n'.format(label))

                ftp_urls = 'fastq_ftp' if 'ebi.ac.uk' in row['fastq_ftp'] else 'submitted_ftp'

                for ftp_path in row[ftp_urls].split(';'):
                    if not ftp_path.startswith('ftp://'):
                        ftp_path = 'ftp://' + ftp_path
                    filename = os.path.basename(ftp_path)
                    assert filename.endswith('.fastq.gz'), 'Expected fastq file {}'.format(ftp_path)

                    if tech_type == 'long':
                        local_name = label + '.part{}.fastq.gz'.format(long_read_parts)
                    else:
                        # rarely, ENA metadata / file reports list three FASTQ files per library
                        # for paired-end reads... example: PRJEB3381
                        mate_num = 0
                        for (mate, ext) in mate1_indicators + mate2_indicators:
                            if filename.endswith(ext):
                                mate_num = mate
                                break
                        if mate_num == 0:
                            _ = logfile.write('WARNING: expecting short, paired-end reads, skipping '
                                              'over file with no mate pair indicator (R1/_1 or R2/_2): {}'.format(ftp_path))
                            continue
                        local_name = label + '_' + row['run_accession'] + '_{}'.format(mate_num) + '.fastq.gz'

                    local_data_folder = os.path.split(request_path)[0]
                    if not requests_created:
                        os.makedirs(request_path, exist_ok=True)
                        requests_created = True

                    local_path = os.path.join(local_data_folder, local_name)
                    request_file = os.path.join(request_path, local_name.replace('.fastq.gz', '.request'))

                    if blacklist and any([hint in ftp_path for hint in blacklist]):
                        _ = logfile.write('File blacklisted: {}\n'.format(local_path))
                        try:
                            os.remove(request_file)
                        except (OSError, IOError) as err:
                            if err.errno != errno.ENOENT:
                                raise RuntimeError('Could not remove blacklisted request file {}: {}'.format(request_file, str(err)))
                        try:
                            os.remove(local_path)
                        except (OSError, IOError) as err:
                            if err.errno != errno.ENOENT:
                                raise RuntimeError('Could not remove blacklisted data file {}: {}'.format(local_path, str(err)))
                        continue

                    _ = logfile.write('Creating request file at {}\n'.format(request_file))
                    try:
                        with open(request_file, 'w') as req:
                            _ = req.write(ftp_path + '\n')
                            _ = req.write(local_path + '\n')
                    except Exception as err:
                        _ = logfile.write('Error: {}'.format(str(err)))
                        raise err
    return requests_created


def find_blacklist_file(given_path):

    absolute_path = None
    search_paths = [
        given_path,
        os.path.join(workflow.basedir, given_path),
        os.path.join(os.getcwd(), given_path)
    ]
    searched_paths = []
    for sp in search_paths:
        if os.path.isfile(sp):
            absolute_path = os.path.abspath(sp)
            break
        searched_paths.append(sp)

    if absolute_path is None:
        raise RuntimeError('Cannot locate blacklist file ({})'.format(searched_paths))

    return absolute_path


checkpoint create_input_data_download_requests:
    input:
        json_dump = rules.master_scrape_data_sources.input,
        tsv_metadata = rules.master_query_repo_sources.input
    output:
        directory('input/{subfolder}/{readset}/requests')
    log:
        'log/input/{subfolder}/{readset}.requests.log'
    params:
        blacklist = lambda wildcards: None if 'file_download_blacklist' not in config else config['file_download_blacklist']
    wildcard_constraints:
        subfolder = '(fastq|bam)'
    run:
        import sys

        try:
            json_dump_files = list(input.json_dump)
        except AttributeError:
            json_dump_files = []

        try:
            tsv_metadata_files = list(input.tsv_metadata)
        except AttributeError:
            tsv_metadata_files = []

        with open(log[0], 'w') as logfile:
            if params.blacklist is None:
                blacklisted_files = set()
            else:
                blacklist_file_path = find_blacklist_file(params.blacklist)
                with open(blacklist_file_path, 'r') as blacklist:
                    blacklisted_files = set(blacklist.read().strip().split())
                _ = logfile.write('Loaded {} files from blacklist annotation\n'.format(len(blacklisted_files)))
            _ = logfile.write('Processing {} JSON dumps\n'.format(len(json_dump_files)))
            json_triggered = create_request_files_from_json(
                                json_dump_files,
                                output[0],
                                blacklisted_files,
                                logfile
                            )
            _ = logfile.write('JSON dump matched request: {}\n'.format(json_triggered))

            _ = logfile.write('Processing {} TSV files\n'.format(len(tsv_metadata_files)))
            tsv_triggered = create_request_files_from_tsv(
                                tsv_metadata_files,
                                output[0],
                                blacklisted_files,
                                logfile
                            )
            _ = logfile.write('TSV file matched request: {}\n'.format(tsv_triggered))

            if json_triggered and tsv_triggered:
                msg = '\nWARNING: request path {} was active for both a JSON dump and a TSV metadata file\n'.format(output[0])
                sys.stderr.write(msg)
                _ = logfile.write(msg)
            elif not json_triggered and not tsv_triggered:
                msg = 'Error: no matching data source for input data request: {}\n'.format(output[0])
                msg += 'Did you forget to load the data source configuration file(s)?\n'
                raise ValueError(msg)
            else:
                pass


rule handle_strandseq_download_requests:
    input:
        'input/fastq/{sseq_reads}/requests/{req_file}.request'
    output:
        'input/fastq/{sseq_reads}/{req_file}.fastq.gz'
    log:
        'log/input/fastq/{sseq_reads}/{req_file}.download.log'
    benchmark:
        'run/input/fastq/{sseq_reads}/{req_file}.download.rsrc'
    wildcard_constraints:
        sseq_reads = CONSTRAINT_STRANDSEQ_SAMPLES
    conda:
         '../environment/conda/conda_shelltools.yml'
    threads: 2
    params:
        script_exec = lambda wildcards: find_script_path('downloader.py', 'utilities'),
        force_copy = lambda wildcards: '--force-local-copy' if bool(config.get('force_local_copy', False)) else ''
    shell:
        '{params.script_exec} --debug {params.force_copy} '
        '--request-file {input} --output {output} '
        '--parallel-conn 1 &> {log}'


rule handle_partial_fastq_download_request:
    input:
        'input/fastq/{readset}/requests/{readset}.{partnum}.request'
    output:
        'input/fastq/{readset}/{readset}.{partnum}.fastq.gz'
    log:
        'log/input/fastq/{readset}/{readset}.{partnum}.download.log'
    benchmark:
        'run/input/fastq/{readset}/{readset}.{partnum}.download.rsrc'
    wildcard_constraints:
        readset = CONSTRAINT_PARTS_FASTQ_INPUT_SAMPLES
    conda:
         '../environment/conda/conda_shelltools.yml'
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    params:
        script_exec = lambda wildcards: find_script_path('downloader.py', 'utilities'),
        force_copy = lambda wildcards: '--force-local-copy' if bool(config.get('force_local_copy', False)) else '',
        parallel_conn = max(1, config['num_cpu_low'])
    shell:
         '{params.script_exec} --debug {params.force_copy}  '
         '--request-file {input} --output {output} '
         '--parallel-conn {params.parallel_conn} &> {log}'


rule handle_partial_pacbio_download_request:
    input:
        'input/bam/{readset}/requests/{readset}.{partnum}.request'
    output:
        'input/bam/{readset}/{readset}.{partnum}.pbn.bam'
    log:
        'log/input/bam/{readset}/{readset}.{partnum}.download.log'
    benchmark:
        'run/input/bam/{readset}/{readset}.{partnum}.download.rsrc'
    wildcard_constraints:
        readset = CONSTRAINT_PARTS_PBN_INPUT_SAMPLES
    conda:
         '../environment/conda/conda_shelltools.yml'
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    params:
        script_exec = lambda wildcards: find_script_path('downloader.py', 'utilities'),
        force_copy = lambda wildcards: '--force-local-copy' if bool(config.get('force_local_copy', False)) else '',
        parallel_conn = max(1, config['num_cpu_low'])
    shell:
         '{params.script_exec} --debug {params.force_copy}  '
         '--request-file {input} --output {output} '
         '--parallel-conn {params.parallel_conn} &> {log}'


def complete_fastq_samples_mock_merger(wildcards):
    """
    This mock-like function exists because Snakemake
    seems to fail recognizing a checkpoint as rule
    dependency if there is no "aggregate" input funtion that
    explicitly calls a "get" on the checkpoint output.
    So, here it is...
    """
    import os
    subfolder = 'fastq'
    readset = wildcards.readset

    requested_input = checkpoints.create_input_data_download_requests.get(subfolder=subfolder, readset=readset).output[0]

    req_file_path = os.path.join(requested_input, readset + '_1000.request')

    return req_file_path


rule handle_single_fastq_download_request:
    input:
        complete_fastq_samples_mock_merger
        #'input/fastq/{readset}/requests/{readset}.request'
    output:
        'input/fastq/{readset}_1000.fastq.gz'
    wildcard_constraints:
        readset = CONSTRAINT_COMPLETE_FASTQ_INPUT_SAMPLES
    log:
        'log/input/fastq/{readset}.download.log'
    benchmark:
        'run/input/fastq/{readset}.download.rsrc'
    conda:
         '../environment/conda/conda_shelltools.yml'
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    params:
        script_exec = lambda wildcards: find_script_path('downloader.py', 'utilities'),
        force_copy = lambda wildcards: '--force-local-copy' if bool(config.get('force_local_copy', False)) else '',
        parallel_conn = max(1, config['num_cpu_low'])
    shell:
         '{params.script_exec} --debug {params.force_copy} '
         '--request-file {input} --output {output} '
         '--parallel-conn {params.parallel_conn} &> {log}'


def select_short_read_request_file(wildcards):
    import os
    subfolder = 'fastq'
    readset = wildcards.readset

    assert wildcards.readset.startswith(wildcards.lib_prefix),\
        'ERROR: sample/readset mismatch for short reads: {} / {}'.format(wildcards.readset, wildcards.lib_prefix)

    requested_input = checkpoints.create_input_data_download_requests.get(subfolder=subfolder, readset=readset).output[0]

    req_file_path = os.path.join(
        requested_input,
        '{}_{}_{}.request'.format(wildcards.lib_prefix, wildcards.lib_id, wildcards.mate))

    return req_file_path


rule handle_short_read_download_request:
    input:
        select_short_read_request_file
    output:
        'input/fastq/{readset}/{lib_prefix}_{lib_id}_{mate}.fastq.gz'
    wildcard_constraints:
        readset = CONSTRAINT_SHORT_READ_INPUT_SAMPLES
    log:
        'log/input/fastq/{readset}/{lib_prefix}_{lib_id}_{mate}.download.log'
    benchmark:
        'run/input/fastq/{readset}/{lib_prefix}_{lib_id}_{mate}.download.rsrc'
    conda:
        '../environment/conda/conda_shelltools.yml'
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    params:
        script_exec = lambda wildcards: find_script_path('downloader.py', 'utilities'),
        force_copy = lambda wildcards: '--force-local-copy' if bool(config.get('force_local_copy', False)) else '',
        parallel_conn = max(1, config['num_cpu_low'])
    shell:
        '{params.script_exec} --debug {params.force_copy} '
        '--request-file {input} --output {output} '
        '--parallel-conn {params.parallel_conn} &> {log}'


def get_bioproject_sample_annotator(bioproject):
    known_projects = {
        'PRJEB12849': sample_annotator_prjeb12849,
        'PRJEB14185': sample_annotator_prjeb14185,
        'PRJNA540705': sample_annotator_prjna540705,
        'PRJEB9396': sample_annotator_prjeb9396,
        'PRJEB36890': sample_annotator_prjeb36890_prjeb31736,
        'PRJEB31736': sample_annotator_prjeb36890_prjeb31736,
    }
    return known_projects[bioproject]


def sample_annotator_prjeb12849(sample_info, individual):

    project = '1kg'
    vendor = 'il'
    assert 'illumina' in sample_info['instrument_platform'].lower()
    model = '25k'
    assert 'hiseq 2500' in sample_info['instrument_model'].lower()
    libname = sample_info['library_name']
    lib_individual, plate, libnum, nuc = libname.split('_')
    assert lib_individual == individual, 'Individual mismatch: {} / {}'.format(libname, individual)
    plate = 'P' + plate.rjust(3, '0')
    assert len(libnum) == 3, 'Unexpected library number: {} / {}'.format(libname, libnum)
    libnum = 'L' + libnum
    sample_id = plate + libnum
    if '_mono' in sample_info['library_name'].lower():
        #read_info = '75pe'
        read_info = '075mo'
    elif '_di' in sample_info['library_name'].lower():
        #read_info = '150pe'
        read_info = '150di'
    else:
        raise ValueError('Unexpected library name: {}'.format(sample_info))
    assert 'paired' in sample_info['library_layout'].lower()
    if individual in sample_info['sample_alias']:
        sample_label = '{}_{}_{}{}-{}_{}'.format(individual, project, vendor, model, read_info, sample_id)
    else:
        sample_label = None
    return sample_label, 'short'


def sample_annotator_prjeb14185(sample_info, individual):
    """
    :param sample_info:
    :param individual:
    :return:
    """
    if sample_info['library_layout'].lower() != 'paired':
        sample_label = None
    else:
        project = 'eriba'
        vendor = 'il'
        read_info = '100pe'
        assert 'illumina' in sample_info['instrument_platform'].lower()
        model = '25k'
        assert 'hiseq 2500' in sample_info['instrument_model'].lower()
        libname = sample_info['library_name']
        libname = libname.replace('_', '-')
        if individual in sample_info['sample_alias']:
            sample_label = '{}_{}_{}{}-{}_{}'.format(individual, project, vendor, model, read_info, libname)
        else:
            sample_label = None
    return sample_label, 'short'


def sample_annotator_prjna540705(sample_info, individual):
    """
    :param sample_info:
    :param individual:
    :return:
    """
    project = 'giab'
    vendor = 'pb'
    read_info = 'ccs'
    assert 'pacbio' in sample_info['instrument_platform'].lower()
    model = 'sq2'
    assert 'sequel ii' in sample_info['instrument_model'].lower()
    if individual in sample_info['sample_alias']:
        sample_label = '{}_{}_{}{}-{}'.format(individual, project, vendor, model, read_info)
    else:
        sample_label = None
    return sample_label, 'long'


def sample_annotator_prjeb9396(sample_info, individual):
    """
    :param sample_info:
    :param individual:
    :return:
    """
    project = '1kg'
    vendor = 'il'
    read_info = '125pe'
    assert 'illumina' in sample_info['instrument_platform'].lower()
    model = '25k'
    assert 'hiseq 2500' in sample_info['instrument_model'].lower()
    assert 'paired' in sample_info['library_layout'].lower()
    if individual in sample_info['sample_alias']:
        sample_label = '{}_{}_{}{}-{}'.format(individual, project, vendor, model, read_info)
    else:
        sample_label = None
    if 'NoIndex' not in sample_info['submitted_ftp']:
        sample_label = None
    return sample_label, 'short'


def sample_annotator_prjeb36890_prjeb31736(sample_info, individual):
    """
    Annotator for NYGC Illumina short read data 698 and 2504 cohort
    Note: the metadata downloaded from ENA contain several sets of FASTQs;
    to be compatible with other approaches, select only those FASTQs dumped
    from submitted CRAM files

    :param sample_info:
    :param individual:
    :return:
    """
    project = '1kg'
    vendor = 'il'
    read_info = '150pe'
    model = 'nvs'
    if individual in sample_info['sample_alias']:
        sample_label = '{}_{}_{}{}-{}'.format(individual, project, vendor, model, read_info)
    else:
        sample_label = None
    if not sample_info['submitted_ftp'].lower().endswith('.cram'):
        sample_label = None
    assert 'illumina' in sample_info['instrument_platform'].lower()
    assert 'novaseq 6000' in sample_info['instrument_model'].lower()
    assert 'paired' in sample_info['library_layout'].lower()
    return sample_label, 'short'