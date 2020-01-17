
localrules: master_handle_data_download


rule master_handle_data_download:
    input:
        rules.master_link_data_sources.output,


checkpoint create_input_data_download_requests:
    input:
        rules.master_scrape_data_sources.input
    output:
        directory('input/{subfolder}/requests')
    wildcard_constraints:
        subfolder = '(fastq|bam)/(complete|partial)[a-z/]*'
    run:
        import os
        import json as json
        import sys as sys
        import collections as col

        all_data_sources = list(input)
        try:
            manual_annotation = config['data_sources']
            assert os.path.isfile(manual_annotation), 'No data sources found at path: {}'.format(manual_annotation)
            all_data_sources.append(manual_annotation)
        except KeyError:
            sys.stderr.write('\nWarning: no manually annotated data sources identified in config\n')

        complete_sources = dict()
        check_key_dups = col.Counter()

        for annotation in all_data_sources:
            with open(annotation, 'r') as paths:
                obj = json.load(paths)
                check_key_dups.update(list(obj.keys()))
                complete_sources.update(obj)

        count_keys = check_key_dups.most_common()
        if count_keys[0][1] > 1:
            raise ValueError('Duplicate keys in data sources annotation: {}'.format(count_keys[0]))

        os.makedirs(output[0], exist_ok=True)
        prefix = os.path.split(output[0])[0]
        os.makedirs(prefix, exist_ok=True)

        for key, file_infos in complete_sources.items():
            if key.startswith(prefix):
                file_prefix = os.path.split(key)[1]
                req_file_path = os.path.join(output[0], file_prefix + '.request')
                try:
                    md5 = file_infos['md5']
                except KeyError:
                    md5 = 'no_md5'
                if os.path.isfile(req_file_path):
                    continue
                with open(req_file_path, 'w') as req_file:
                    _ = req_file.write(file_infos['remote_path'] + '\n')
                    _ = req_file.write(file_infos['local_path'] + '\n')
                    _ = req_file.write(md5 + '\n')
    # end of rule


def select_strandseq_bioproject_source(sts_reads):
    """
    :param sts_reads:
    :return:
    """
    individual = sts_reads.split('_')[0]
    sample_desc = config['sample_description_' + individual]
    data_sources = sample_desc['data_sources']

    bioproject = None
    for record in data_sources:
        if 'strandseq' not in record:
            continue
        record_spec = record['strandseq']
        if record_spec['readset'].startswith(sts_reads):
            # not sure if prefix matching necessary... should always be full match
            if 'bioproject' not in record_spec:
                raise ValueError('No bioproject annotated in config for Strand-seq readset: {}'.format(sts_reads))
            bioproject = record_spec['bioproject']
    if bioproject is None:
        raise ValueError('No matching bioproject annotated for Strand-seq data: {}'.format(sts_reads))
    return bioproject


rule download_bioproject_metadata:
    output:
        'input/bioprojects/{sts_reads}.metadata.tsv'
    log:
       'log/input/bioprojects/{sts_reads}.metadata.log'
    conda:
         '../environment/conda/conda_shelltools.yml'
    params:
        script_dir = config['script_dir'],
        accession_number = lambda wildcards: select_strandseq_bioproject_source(wildcards.sts_reads)
    shell:
        '{params.script_dir}/utilities/downloader.py --debug '
        '--ena-file-report '
        ' "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession='
        '{params.accession_number}'
        '&result=read_run&fields='
        'study_accession,sample_accession,secondary_sample_accession'
        ',experiment_accession,run_accession,submission_accession'
        ',tax_id,scientific_name,instrument_platform,instrument_model'
        ',library_name,library_layout,library_strategy,library_source'
        ',library_selection,read_count,center_name,study_title,fastq_md5'
        ',study_alias,experiment_alias,run_alias,fastq_ftp,submitted_ftp'
        ',sample_alias,sample_title'
        '\&format=tsv\&download=txt" '
        '--output {output} &> {log}'


checkpoint create_bioproject_download_requests:
    input:
        'input/bioprojects/{sts_reads}.metadata.tsv'
    output:
        directory('input/fastq/strand-seq/{sts_reads}/requests')
    log:
        'log/input/fastq/strand-seq/{sts_reads}.requests.log'
    run:
        import os
        import csv

        bioproject = select_strandseq_bioproject_source(wildcards.sts_reads)
        individual = wildcards.sts_reads.split('_')[0]

        with open(log[0], 'w') as logfile:
            annotator = get_bioproject_sample_annotator(bioproject)
            with open(input[0], 'r') as table:
                rows = csv.DictReader(table, delimiter='\t')
                for row in rows:
                    if individual not in row['sample_alias']:
                        continue
                    label = annotator(row, individual)
                    if label is None:
                        continue
                    _ = logfile.write('Processing data for sample {}\n'.format(label))

                    for ftp_path in row['fastq_ftp'].split(';'):
                        if not ftp_path.startswith('ftp://'):
                            ftp_path = 'ftp://' + ftp_path
                        filename = os.path.basename(ftp_path)
                        local_name = label + '_' + filename

                        local_fastq_folder = os.path.split(output[0])[0]
                        local_req_folder = output[0]
                        os.makedirs(local_req_folder, exist_ok=True)

                        local_fastq_path = os.path.join(local_fastq_folder, local_name)
                        request_file = os.path.join(local_req_folder, local_name.replace('.fastq.gz', '.request'))

                        _ = logfile.write('Creating request file at {}\n'.format(request_file))
                        try:
                            with open(request_file, 'w') as req:
                                _ = req.write(ftp_path + '\n')
                                _ = req.write(local_fastq_path + '\n')
                        except Exception as err:
                            _ = logfile.write('Error: {}'.format(str(err)))
                            raise err
            _ = logfile.write('Done - create_bioproject_download_requests')
    # end of rule


rule handle_strandseq_download_requests:
    input:
        'input/fastq/strand-seq/{sts_reads}/requests/{sample}_{run_id}.request'
    output:
        'input/fastq/strand-seq/{sts_reads}/{sample}_{run_id}.fastq.gz'
    log:
        'log/input/fastq/strand-seq/{sts_reads}/{sample}_{run_id}.download.log'
    benchmark:
        'run/input/fastq/strand-seq/{sts_reads}/{sample}_{run_id}.download.rsrc'
    conda:
         '../environment/conda/conda_shelltools.yml'
    threads: 2
    params:
        script_dir = config['script_dir']
    shell:
        '{params.script_dir}/utilities/downloader.py --debug '
        '--request-file {input} --output {output} '
        '--parallel-conn 1 &> {log}'


rule handle_partial_fastq_download_request:
    input:
        'input/fastq/partial/{split_type}/requests/{req_sample}.{partnum}.request'
    output:
        'input/fastq/partial/{split_type}/{req_sample}.{partnum}.fastq.gz'
    log:
        'log/input/fastq/partial/{split_type}/{req_sample}.{partnum}.download.log'
    benchmark:
        'run/input/fastq/partial/{split_type}/{req_sample}.{partnum}.download.rsrc'
    wildcard_constraints:
        split_type = '(parts|chunks)',
        req_sample = CONSTRAINT_PARTS_FASTQ_INPUT_SAMPLES
    conda:
         '../environment/conda/conda_shelltools.yml'
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    params:
        script_dir = config['script_dir'],
        parallel_conn = config['num_cpu_low'] - 1
    shell:
         '{params.script_dir}/utilities/downloader.py --debug '
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
    subfolder = 'fastq/complete'

    requested_input = checkpoints.create_input_data_download_requests.get(subfolder=subfolder).output[0]

    base_path = os.path.join('input', subfolder)
    request_path = os.path.join(base_path, 'requests')

    sample = wildcards.req_sample

    req_file_path = os.path.join(request_path, sample + '.request')
    assert os.path.isfile(req_file_path), 'Path not file after checkpoint execution: {}'.format(req_file_path)

    return req_file_path


rule handle_complete_fastq_download_request:
    input:
        complete_fastq_samples_mock_merger
        #'input/fastq/complete/requests/{sample}.request'
    output:
        'input/fastq/complete/{req_sample}_1000.fastq.gz'
    wildcard_constraints:
        req_sample = CONSTRAINT_COMPLETE_FASTQ_INPUT_SAMPLES
    log:
        'log/input/fastq/complete/{req_sample}.download.log'
    benchmark:
        'run/input/fastq/complete/{req_sample}.download.rsrc'
    conda:
         '../environment/conda/conda_shelltools.yml'
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    params:
        script_dir = config['script_dir'],
        parallel_conn = config['num_cpu_low'] - 1
    shell:
         '{params.script_dir}/utilities/downloader.py --debug '
         '--request-file {input} --output {output} '
         '--parallel-conn {params.parallel_conn} &> {log}'


rule handle_partial_pbn_bam_download_request:
    input:
        'input/bam/partial/{split_type}/requests/{req_sample}.{partnum}.request'
    output:
        'input/bam/partial/{split_type}/{req_sample}.{partnum}.pbn.bam'
    log:
        'log/input/bam/partial/{split_type}/{req_sample}.{partnum}.download.log'
    benchmark:
        'run/input/bam/partial/{split_type}/{req_sample}.{partnum}.download.rsrc'
    conda:
        '../environment/conda/conda_shelltools.yml'
    wildcard_constraints:
        split_type = '(parts|chunks)',
        req_sample = CONSTRAINT_PARTS_PBN_INPUT_SAMPLES
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: 8 * attempt
    params:
        script_dir = config['script_dir'],
        parallel_conn = config['num_cpu_low'] - 1
    shell:
         '{params.script_dir}/utilities/downloader.py --debug '
         '--request-file {input} --output {output} '
         '--parallel-conn {params.parallel_conn} &> {log}'


def get_bioproject_sample_annotator(bioproject):
    known_projects = {'PRJEB12849': sample_annotator_PRJEB12849}
    return known_projects[bioproject]


def sample_annotator_PRJEB12849(sample_info, individual):

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
    return sample_label