
localrules: master_scrape_data_sources


def read_configured_data_sources():
    """
    :return:
    """

    COMMAND_SCAN_PATH_CALL = '{{params.script_exec}} --debug --output {{output}} '
    COMMAND_SCAN_PATH_CALL += '--server {server} --data-source {data_source} '
    COMMAND_SCAN_PATH_CALL += '--collect-files {collect_files} --sort-into {sort_into} '
    COMMAND_SCAN_PATH_CALL += '--file-infix {file_infix} {file_suffix} '
    COMMAND_SCAN_PATH_CALL += '{local_path_suffix} {fix_tech} '
    COMMAND_SCAN_PATH_CALL += '{assume_pacbio_native} {assume_clr_subreads} {assume_paired_reads} '
    COMMAND_SCAN_PATH_CALL += '&> {{log}}'

    import os
    import sys

    default_parameters = {
        'file_suffix': '',
        'local_path_suffix': '',
        'fix_tech': '',
        'assume_pacbio_native': '',
        'assume_clr_subreads': '',
        'assume_paired_reads': ''
    }

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
            source_syscalls.append(source_call)

    if len(source_outputs) != len(set(source_outputs)):
        sys.stderr.write('\nERROR: at least two data sources specify the same output file: {}\n'.format(sorted(source_outputs)))
        raise ValueError('Duplicate data source output file.')

    if not source_names:
        sys.stderr.write('\nWARNING: no data sources configured\n')

    return source_names, source_outputs, source_syscalls


DATA_SOURCE_NAMES, DATA_SOURCE_OUTPUTS, DATA_SOURCE_CALLS = read_configured_data_sources()


rule master_scrape_data_sources:
    input:
        DATA_SOURCE_OUTPUTS


for name, outfile, syscall in zip(DATA_SOURCE_NAMES, DATA_SOURCE_OUTPUTS, DATA_SOURCE_CALLS):
    rule:
        output:
            outfile
        log:
            os.path.join('log', outfile.replace('.json', '.log'))
        message: 'Processing data source: {}'.format(name)
        conda:
            '../environment/conda/conda_pyscript.yml'
        params:
            script_exec = lambda wildcards: find_script_path('scan_remote_path.py')
        shell:
            syscall
