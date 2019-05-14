
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
                remote_url = v['url']
                local_name = v['name']
                if reftype == 'genome':
                    use_genome = config['use_ref_genome']
                    if not local_name.startswith(use_genome):
                        continue
                if 'path' in v:
                    local_path = os.path.join(v['path'], local_name)
                else:
                    local_path = os.path.join('references', local_name)
                rule:
                    output:
                        local_path
                    log: 'log/references/DL_{}.log'.format(local_name.rsplit('.', 1)[0])
                    threads: 4
                    message: 'Downloading reference file {} / {}'.format(reftype, local_name)
                    run:
                        exec = 'aria2c -s {threads} -x {threads}'
                        exec += ' -o {output}'
                        exec += ' {}'.format(remote_url)
                        exec += ' &> {log}'
                        shell(exec)