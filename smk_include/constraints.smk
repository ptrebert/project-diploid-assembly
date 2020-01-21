
def build_input_data_wildcard_constraint(input_type, readset_selectors):
    """
    :param input_type:
    :param readset_selectors:
    :return:
    """
    import sys
    selected_readsets = []

    for config_key, settings in config.items():
        if not config_key.startswith('sample_description'):
            continue
        try:
            sources = settings['data_sources']
        except KeyError:
            sys.stderr.write('\nWARNING: no data sources configured for sample: {}\n'.format(settings['individual']))
            continue
        for source in sources:
            if input_type not in source:
                continue
            source_spec = source[input_type]
            if 'readset' not in source_spec:
                sys.stderr.write('\nWARNING: no key "readset" found in data '
                                 'source specification: {} / {} / {}\n'.format(settings['individual'], input_type, readset_selectors))
                continue
            select = True
            for select_key, select_values in readset_selectors.items():
                assert isinstance(select_values, list), \
                    'Select values must be given as a list of possible values: {}'.format(readset_selectors)
                if select_key not in source_spec:
                    select = False
                    # missing key means this readset does not
                    # have the desired attribute
                    break
                if source_spec[select_key] not in select_values:
                    # note that values should always be a
                    # list to avoid matching substrings like that
                    select = False
                    break
            if select:
                selected_readsets.append(source_spec['readset'])

    if not selected_readsets:
        sys.stderr.write('\nWARNING: no wildcard constraint created '
                         'for readset selectors: {} / {}\n'.format(input_type, readset_selectors))
        wildcard_regexp = '^$'
    else:
        wildcard_regexp = '(' + '|'.join(sorted(selected_readsets)) + ')'

    return wildcard_regexp


CONSTRAINT_ALL_PBN_INPUT_SAMPLES = build_input_data_wildcard_constraint(
    input_type='long_reads',
    readset_selectors={
        'data_type': ['pacbio_native'],
        'load_type': ['parts', 'complete']
    }
)

CONSTRAINT_COMPLETE_PBN_INPUT_SAMPLES = build_input_data_wildcard_constraint(
    input_type='long_reads',
    readset_selectors={
        'data_type': ['pacbio_native'],
        'load_type': ['complete']
    }
)

CONSTRAINT_PARTS_PBN_INPUT_SAMPLES = build_input_data_wildcard_constraint(
    input_type='long_reads',
    readset_selectors={
        'data_type': ['pacbio_native'],
        'load_type': ['parts']
    }
)

CONSTRAINT_ALL_FASTQ_INPUT_SAMPLES = build_input_data_wildcard_constraint(
    input_type='long_reads',
    readset_selectors={
        'data_type': ['fastq'],
        'load_type': ['parts', 'complete']
    }
)

CONSTRAINT_COMPLETE_FASTQ_INPUT_SAMPLES = build_input_data_wildcard_constraint(
    input_type='long_reads',
    readset_selectors={
        'data_type': ['fastq'],
        'load_type': ['complete']
    }
)

CONSTRAINT_PARTS_FASTQ_INPUT_SAMPLES = build_input_data_wildcard_constraint(
    input_type='long_reads',
    readset_selectors={
        'data_type': ['fastq'],
        'load_type': ['parts']
    }
)