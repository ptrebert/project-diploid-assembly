
localrules: master_link_data_sources

rule master_link_data_sources:
    """
    This is the place to inject external data
    into the pipeline via symlinking
    """
    output: