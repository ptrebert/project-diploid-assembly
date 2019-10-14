
CMD_DL_COMPRESSED_PARALLEL = 'aria2c --out={{output}} --file-allocation=none -s 4 -x 4 {remote_path} &>> {{log}}'

CMD_DL_COMPRESSED_SINGLE = 'wget -O {{output}} {remote_path} &>> {{log}}'

CMD_DL_UNCOMPRESSED_SINGLE = 'wget --quiet -O /dev/stdout {remote_path} 2>> {{log}} | gzip > {{output}}'

localrules: master_handle_reference_download, create_reference_download_request

rule master_handle_reference_download:
    input:


rule create_reference_download_request:
    output:
        'references/{subfolder}/{reference}.request'
    run:
        import json as json

        ref_sources_path = config['reference_sources']
        assert os.path.isfile(ref_sources_path), 'No reference sources configured'

        with open(ref_sources_path, 'r') as annotation:
            sources = json.load(annotation)

        file_info = sources[wildcards.reference]

        try:
            md5 = file_info['md5']
        except KeyError:
            md5 = 'no_md5'

        with open(output[0], 'w') as req_file:
            _ = req_file.write(file_info['remote_path'] + '\n')
            _ = req_file.write(file_info['local_path'] + '\n')
            _ = req_file.write(md5 + '\n')
    # end of rule


rule handle_raw_fasta_reference_download_request:
    input:
        'references/downloads/{reference}.request'
    output:
        'references/downloads/{reference}.fa.gz'
    log:
        'log/references/downloads/{reference}.download.log'
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

rule handle_gff_reference_download_request:
    input:
        'references/downloads/{reference}.request'
    output:
        'references/downloads/{reference}.gff3.gz'
    log:
        'log/references/downloads/{reference}.download.log'
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