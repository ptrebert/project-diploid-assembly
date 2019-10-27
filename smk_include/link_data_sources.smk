
localrules: master_link_data_sources

rule master_link_data_sources:
    """
    This is the place to inject external data
    into the pipeline via symlinking
    Attention: Snakemake may delete "output" files upon a re-run,
    so, presumably, they should be "protected" AND "ancient"
    to (hopefully) ensure that the pipeline runs smoothly
    """
    output: