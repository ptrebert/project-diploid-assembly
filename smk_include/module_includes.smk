
# modules w/o other dependencies
include: 'constraints.smk'
include: 'aux_utilities.smk'
include: 'environments.smk'
include: 'link_data_sources.smk'
include: 'scrape_data_sources.smk'
include: 'query_data_repos.smk'
include: 'handle_reference_download.smk'

# input preparation stage, one or two dependencies to above modules
include: 'handle_data_download.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'

# actual pipeline processing steps
include: 'variant_calling.smk'
include: 'integrative_phasing.smk'

include: 'strandseq_dga_split.smk'
include: 'strandseq_dga_joint.smk'

include: 'collect_statistics.smk'
include: 'run_alignments.smk'
include: 'run_assemblies.smk'
include: 'prepare_custom_references.smk'

include: 'run_polishing.smk'

include: 'haploid_assembly_clustering.smk'

include: 'create_plots.smk'
include: 'eval_known_reference.smk'

include: 'targets.smk'