
CMD_DL_COMPRESSED_PARALLEL = 'aria2c --out={{output}} --file-allocation=none -s 4 -x 4 {remote_path} &>> {{log}}'

CMD_DL_COMPRESSED_SINGLE = 'wget -O {{output}} {remote_path} &>> {{log}}'

CMD_DL_UNCOMPRESSED_SINGLE = 'wget --quiet -O /dev/stdout {remote_path} 2>> {{log}} | gzip > {{output}}'

localrules: master_handle_data_download, create_input_data_download_request

rule master_handle_data_download:
    input:
        expand('input/fastq/partial/parts/{sample}.part{partnum}.fastq.gz',
                sample=['HG00733_hpg_ontpm-ul'],
                partnum=[1, 2, 3])


rule create_input_data_download_request:
    input:
        config['data_sources']
    output:
        'input/{subfolder}/{sample}.request'
    run:
        import json as json

        with open(input[0], 'r') as annotation:
            sources = json.load(annotation)

        file_info = sources[wildcards.sample]

        try:
            md5 = file_info['md5']
        except KeyError:
            md5 = 'no_md5'

        with open(output[0], 'w') as req_file:
            _ = req_file.write(file_info['remote_path'] + '\n')
            _ = req_file.write(file_info['local_path'] + '\n')
            _ = req_file.write(md5 + '\n')
    # end of rule


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

        with open(log[0], 'w') as logfile:
            _ = logfile.write('Handling download request' + '\n')
            _ = logfile.write('CMD: {}'.format(exec) + '\n\n')
            if local_path != output[0]:
                _ = logfile.write('ERROR - output mismatch: {} vs {}'.format(local_path, output[0]))
                raise RuntimeError

        shell(exec)
    # end of rule