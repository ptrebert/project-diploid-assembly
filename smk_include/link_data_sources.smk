
localrules: master_link_data_sources

rule master_link_data_sources:
    """
    This is the place to inject external data
    into the pipeline via symlinking
    """
    input:
        ancient(config.get('link_data_input', []))
    output:
        config.get('link_data_output', [])
    run:
        input_files = list(input)
        output_links = list(output)

        if len(input_files) != len(output_links):
            raise RuntimeError('Cannot inject data via sym linking, no 1-to-1 correspondence '
                               'between input and output: {} vs {}'.format(len(input_files), len(output_links)))

        import os
        for input_file, output_link in zip(input_files, output_links):
            assert os.path.isfile(input_file), 'Invalid path to input file for linking: {}'.format(input_file)
            os.makedirs(os.path.dirname(output_link), exist_ok=True)
            os.symlink(input_file, output_link)