
CMD_DL_COMPRESSED_PARALLEL = 'aria2c --out={{output}} --file-allocation=none -s 4 -x 4 {remote_path} &>> {{log}}'

CMD_DL_COMPRESSED_SINGLE = 'wget -O {{output}} {remote_path} &>> {{log}}'

CMD_DL_UNCOMPRESSED_SINGLE = 'wget --quiet -O /dev/stdout {remote_path} 2>> {{log}} | gzip > {{output}}'

localrules: master_handle_data_download, \
            create_input_data_download_request, \
            download_bioproject_metadata, \
            create_bioproject_download_requests

rule master_handle_data_download:
    input:
        expand('input/fastq/partial/parts/{sample}.part{partnum}.fastq.gz',
                sample=['HG00733_hpg_ontpm-ul'],
                partnum=[1, 2, 3]),

        expand('input/fastq/complete/{sample}_1000.fastq.gz',
                sample=['HG00733_sra_pbsq1-clr']),
        'input/bioprojects/PRJEB12849.tsv',
        'input/fastq/strand-seq/HG00733_PRJEB12849/requests'


rule create_input_data_download_request:
    output:
        'input/{subfolder}/{sample}.request'
    wildcard_constraints:
        subfolder = '^(fastq|bam).*'
    run:
        import json as json

        data_sources_path = config['data_sources']
        assert os.path.isfile(data_sources_path), 'No data sources configured'

        with open(data_sources_path, 'r') as annotation:
            sources = json.load(annotation)

        file_info = sources[wildcards.sample]

        try:
            md5 = file_info['md5']
        except KeyError:
            md5 = 'no_md5'

        if not os.path.isfile(output[0]):
            with open(output[0], 'w') as req_file:
                _ = req_file.write(file_info['remote_path'] + '\n')
                _ = req_file.write(file_info['local_path'] + '\n')
                _ = req_file.write(md5 + '\n')
    # end of rule


rule download_bioproject_metadata:
    output:
        'input/bioprojects/{bioproject}.tsv'
    run:
        load_url = 'https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession='
        load_url += '{accession}'
        load_url += '&result=read_run&fields='
        load_url += 'study_accession,sample_accession,secondary_sample_accession'
        load_url += ',experiment_accession,run_accession,submission_accession'
        load_url += ',tax_id,scientific_name,instrument_platform,instrument_model'
        load_url += ',library_name,library_layout,library_strategy,library_source'
        load_url += ',library_selection,read_count,center_name,study_title,fastq_md5'
        load_url += ',study_alias,experiment_alias,run_alias,fastq_ftp,submitted_ftp'
        load_url += ',sample_alias,sample_title'
        load_url += '&download=txt'

        tmp = load_url.format(**{'accession': wildcards.bioproject})
        if os.path.isfile(output[0]):
            exec = 'touch {output}'
        else:
            exec = 'wget --quiet -O {{output}} "{}"'.format(tmp)
        shell(exec)


checkpoint create_bioproject_download_requests:
    input:
        'input/bioprojects/{bioproject}.tsv'
    output:
        directory('input/fastq/strand-seq/{individual}_{bioproject}/requests')
    log:
        'log/input/fastq/strand-seq/{individual}_{bioproject}.requests.log'
    run:
        import csv
        #import pdb
        with open(log[0], 'w') as logfile:
            annotator = get_bioproject_sample_annotator(wildcards.bioproject)
            with open(input[0], 'r') as table:
                rows = csv.DictReader(table, delimiter='\t')
                for row in rows:
                    label = annotator(row, wildcards.individual)
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

                        if request_file.endswith('_2.request'):
                            # skip as SaarClust does not used mate2 info
                            continue
                        _ = logfile.write('Creating request file at {}\n'.format(request_file))
                        try:
                            with open(request_file, 'w') as req:
                                _ = req.write(ftp_path + '\n')
                                _ = req.write(local_fastq_path + '\n')
                        except Exception as err:
                            _ = logfile.write('Error: {}'.format(str(err)))
                            raise err
            #pdb.set_trace()
            _ = logfile.write('Done - create_bioproject_download_requests')
    # end of rule


rule handle_strandseq_download_requests:
    input:
        'input/fastq/strand-seq/{individual}_{bioproject}/requests/{sample}_{run_id}.request'
    output:
        'input/fastq/strand-seq/{individual}_{bioproject}/{sample}_{run_id}.fastq.gz'
    log:
        'log/input/fastq/strand-seq/{individual}_{sample}_{bioproject}_{run_id}.download.log'
    threads: 2
    run:
        with open(input[0], 'r') as req_file:
            remote_path = req_file.readline().strip()
            local_path = req_file.readline().strip()

        if remote_path.endswith('.gz'):
            exec = CMD_DL_COMPRESSED_PARALLEL.format(**{'remote_path': remote_path})
        else:
            exec = CMD_DL_UNCOMPRESSED_SINGLE.format(**{'remote_path': remote_path})

        if not os.path.isfile(output[0]):
            with open(log[0], 'w') as logfile:
                _ = logfile.write('Handling download request' + '\n')
                _ = logfile.write('CMD: {}'.format(exec) + '\n\n')
                if local_path != output[0]:
                    _ = logfile.write('ERROR - output mismatch: {} vs {}'.format(local_path, output[0]))
                    raise RuntimeError

            shell(exec)
    # end of checkpoint


rule handle_partial_fastq_download_request:
    input:
        'input/fastq/partial/{split_type}/{sample_split}.request'
    output:
        protected('input/fastq/partial/{split_type}/{sample_split}.fastq.gz')
    log:
        'log/input/fastq/partial/{split_type}/{sample_split}.download.log'
    wildcard_constraints:
        split_type = '(parts|chunks)'
    threads: 2  # compromise between wget and aria2c
    run:
        with open(input[0], 'r') as req_file:
            remote_path = req_file.readline().strip()
            local_path = req_file.readline().strip()

        if remote_path.endswith('.gz'):
            exec = CMD_DL_COMPRESSED_PARALLEL.format(**{'remote_path': remote_path})
        else:
            exec = CMD_DL_UNCOMPRESSED_SINGLE.format(**{'remote_path': remote_path})

        if not os.path.isfile(output[0]):
            with open(log[0], 'w') as logfile:
                _ = logfile.write('Handling download request' + '\n')
                _ = logfile.write('CMD: {}'.format(exec) + '\n\n')
                if local_path != output[0]:
                    _ = logfile.write('ERROR - output mismatch: {} vs {}'.format(local_path, output[0]))
                    raise RuntimeError

            shell(exec)
    # end of rule


rule handle_complete_fastq_download_request:
    input:
        'input/fastq/complete/{sample}.request'
    output:
        protected('input/fastq/complete/{sample}_1000.fastq.gz')
    wildcard_constraints:
        sample = '(' + '|'.join(config['complete_samples']) + ')'
    log:
        'log/input/fastq/complete/{sample}.download.log'
    threads: 2  # compromise between wget and aria2c
    run:
        with open(input[0], 'r') as req_file:
            remote_path = req_file.readline().strip()
            local_path = req_file.readline().strip()

        if remote_path.endswith('.gz'):
            exec = CMD_DL_COMPRESSED_PARALLEL.format(**{'remote_path': remote_path})
        else:
            exec = CMD_DL_UNCOMPRESSED_SINGLE.format(**{'remote_path': remote_path})

        if not os.path.isfile(output[0]):
            with open(log[0], 'w') as logfile:
                _ = logfile.write('Handling download request' + '\n')
                _ = logfile.write('CMD: {}'.format(exec) + '\n\n')
                if local_path != output[0]:
                    _ = logfile.write('ERROR - output mismatch: {} vs {}'.format(local_path, output[0]))
                    raise RuntimeError

            shell(exec)
    # end of rule


def get_bioproject_sample_annotator(bioproject):
    known_projects = {'PRJEB12849': sample_annotator_PRJEB12849}
    return known_projects[bioproject]


def sample_annotator_PRJEB12849(sample_info, individual):

    project = '1kg'
    vendor = 'il'
    assert 'illumina' in sample_info['instrument_platform'].lower()
    model = '25k'
    assert 'hiseq 2500' in sample_info['instrument_model'].lower()
    if '_mono' in sample_info['library_name'].lower():
        #read_info = '75pe'
        read_info = 'npe'
    elif '_di' in sample_info['library_name'].lower():
        #read_info = '150pe'
        read_info = 'npe'
    else:
        raise ValueError('Unexpected library name: {}'.format(sample_info))
    assert 'paired' in sample_info['library_layout'].lower()
    if individual in sample_info['sample_alias']:
        sample_label = '{}_{}_{}{}-{}'.format(individual, project, vendor, model, read_info)
    else:
        sample_label = None
    return sample_label