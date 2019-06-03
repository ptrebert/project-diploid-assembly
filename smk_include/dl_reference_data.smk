
# unnamed rule to download reference
# annotation such as genomes
# all of this lands in /references by default

AUTOLOAD_REFERENCE_DATA = ['ref_data', 'ref_variants']

for config_entry in AUTOLOAD_REFERENCE_DATA:
    for reftype, records in config[config_entry].items():
        for record in records:
            for v in record.values():
                if 'skip' in v:
                    if bool(v['skip']):
                        continue
                local_name = v['name']
                if reftype == 'genome':
                    use_genome = config['use_ref_genome']
                    if not local_name.startswith(use_genome):
                        continue
                if 'path' in v:
                    local_path = os.path.join(v['path'], local_name)
                else:
                    local_path = os.path.join('references', local_name)
                use_threads = 4
                rule:
                    output:
                        local_path
                    log: 'log/references/DL_{}.log'.format(local_name.rsplit('.', 1)[0])
                    threads: use_threads
                    message: 'Downloading reference file {} / {}'.format(reftype, local_name)
                    params:
                        remote_url = v['url']
                    run:
                        source_gzip = params.remote_url.endswith('.gz')
                        target_gzip = local_path.endswith('.gz')
                        load_parallel = (source_gzip and target_gzip) or (not source_gzip and not target_gzip)
                        if load_parallel:
                            exec = 'aria2c -s {threads} -x {threads}'
                            exec += ' -o {output}'
                            exec += ' {params.remote_url}'
                            exec += ' &> {log}'
                        elif source_gzip and (not target_gzip):
                            exec = 'wget --quiet -O -'  # write to stdout
                            exec += ' {params.remote_url} 2> {log}'
                            exec += ' | gunzip -c > {output}'
                            exec += ' 2>> {log}'
                        elif (not source_gzip) and target_gzip:
                            exec = 'wget --quiet -O -'  # write to stdout
                            exec += ' {params.remote_url} 2> {log}'
                            exec += ' | gzip -c > {output}'
                            exec += ' 2>> {log}'
                        else:
                            raise ValueError('Cannot handle combo of remote'
                                             ' and local path: {} / {}'.format(params.remote_url, local_path))
                        shell(exec)