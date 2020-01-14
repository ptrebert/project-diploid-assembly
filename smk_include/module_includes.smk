
# modules w/o other dependencies
include: 'constraints.smk'
include: 'environments.smk'
include: 'aux_utilities.smk'
include: 'link_data_sources.smk'
include: 'scrape_data_sources.smk'
include: 'handle_reference_download.smk'

# input preparation stage, one or two dependencies to above modules
include: 'handle_data_download.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'

include: 'variant_calling.smk'
include: 'integrative_phasing.smk'

include: 'strandseq_dga_joint.smk'
include: 'strandseq_dga_split.smk'

# actual pipeline processing steps
include: 'collect_statistics.smk'
include: 'run_alignments.smk'
include: 'run_assemblies.smk'
include: 'prepare_custom_references.smk'

include: 'arrow_polishing.smk'
include: 'racon_polishing.smk'

include: 'eval_known_reference.smk'

include: 'targets.smk'