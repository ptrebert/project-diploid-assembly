
import os
import sys

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
        sys.stderr.write('\nWARNING: no data sources configured\n')

    if not len(source_names) == len(source_outputs) == len(source_syscalls):
        raise RuntimeError('Length mismatch for data sources: \n'
                           ' => {} \n => {} \n => {}'.format(source_names, source_outputs, source_syscalls))

    return source_names, source_outputs, source_syscalls


DATA_SOURCE_NAMES, DATA_SOURCE_OUTPUTS, DATA_SOURCE_CALLS = read_configured_data_sources()

DATA_SOURCE_TO_CALL = dict((os.path.basename(source).rsplit('.', 1)[0], call) for source, call in zip(DATA_SOURCE_OUTPUTS, DATA_SOURCE_CALLS))


rule master_scrape_data_sources:
    input:
        DATA_SOURCE_OUTPUTS


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
