
localrules: master_handle_reference_download

rule master_handle_reference_download:
    input:
        []


rule create_reference_download_request:
    output:
        'references/{subfolder}/{reference}.request'
    run:
        import os
        ref_sources_info = config['reference_data_sources']

        file_info = ref_sources_info[wildcards.reference]

        request_path = os.path.join('references', wildcards.subfolder)
        if not file_info['local_path'].startswith(request_path):
            raise ValueError('Reference file mismatch: '
                             'requested path is {}, '
                             'annotated is {}'.format(request_path, file_info['local_path']))

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


rule handle_gff_reference_download_request:
    input:
        'references/downloads/{reference}.request'
    output:
        'references/downloads/{reference}.gff3.gz'
    log:
        'log/references/downloads/{reference}.download.log'
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