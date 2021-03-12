
import os
import sys

localrules: master_query_repo_sources


def read_configured_repo_sources():
    """
    :return:
    """

    COMMAND_QUERY_REPO_CALL = '{script_exec} --debug --output {output} '
    COMMAND_QUERY_REPO_CALL += '{repo} {data_source} '
    COMMAND_QUERY_REPO_CALL += '&> {log}'

    possible_repos = {
        'ena': '--ena-file-report'
    }

    script_path = find_script_path('downloader.py', 'utilities')

    source_names = []
    source_outputs = []
    source_syscalls = []

    for key, values in config.items():
        if not key.startswith('sample_description'):
            continue
        sample_name = key.split('_', 2)[-1]

        data_sources = values['data_sources']

        for data_source in data_sources:
            for source_attr in data_source.values():
                source_type = source_attr.get('source_type', None)
                if source_type not in possible_repos:
                    continue
                try:
                    bioproject = source_attr['bioproject']
                except KeyError as error:
                    sys.stderr.write('\nERROR: repository data source does not have a '
                                     'bioproject accession annotated: {}\n'.format(str(error)))
                    raise error
                readset = source_attr['readset']
                assert readset.startswith(sample_name), 'Sample name not in ' \
                                                        'readset name: {} / {}'.format(sample_name, readset)

                repo_source_name = '_'.join([source_type, bioproject, readset])
                repo_source_output = os.path.join('input/data_sources/{}.metadata.tsv'.format(readset))
                repo_source_log = os.path.join('log', repo_source_output.replace('.tsv', '.log'))

                repo_source_call = COMMAND_QUERY_REPO_CALL.format(**{
                    'script_exec': script_path,
                    'output': repo_source_output,
                    'repo': possible_repos[source_type],
                    'data_source': bioproject,
                    'log': repo_source_log
                })

                source_names.append(repo_source_name)
                source_outputs.append(repo_source_output)
                source_syscalls.append(repo_source_call)

    if len(source_outputs) != len(set(source_outputs)):
        sys.stderr.write('\nERROR: at least two data sources specify the same output file: {}\n'.format(sorted(source_outputs)))
        raise ValueError('Duplicate data source output file.')

    if not source_names:
        if bool(config.get('show_warnings', False)):
            sys.stderr.write('\nWARNING: no repository data sources configured\n')

    if not len(source_names) == len(source_outputs) == len(source_syscalls):
        raise RuntimeError('Length mismatch for data sources: \n'
                           ' => {} \n => {} \n => {}'.format(source_names, source_outputs, source_syscalls))

    return source_names, source_outputs, source_syscalls


if bool(config.get('use_legacy_data_scraping', False)):

    REPO_SOURCE_NAMES, REPO_SOURCE_OUTPUTS, REPO_SOURCE_CALLS = read_configured_repo_sources()

    REPO_SOURCE_TO_CALL = dict((os.path.basename(source).rsplit('.', 2)[0], call) for source, call in zip(REPO_SOURCE_OUTPUTS, REPO_SOURCE_CALLS))

else:
    REPO_SOURCE_OUTPUTS = list()
    REPO_SOURCE_TO_CALL = dict()


rule master_query_repo_sources:
    input:
        REPO_SOURCE_OUTPUTS


rule query_repo_sources:
    """
    """
    output:
        'input/data_sources/{data_source}.metadata.tsv'
    log:
        'log/input/data_sources/{data_source}.metadata.log'
    message: 'Processing data source: {output}'
    conda:
        '../environment/conda/conda_shelltools.yml'
    params:
        query_call = lambda wildcards: REPO_SOURCE_TO_CALL[wildcards.data_source]
    shell:
        '{params.query_call}'
