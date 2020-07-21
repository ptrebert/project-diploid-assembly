
include: 'prep_custom_references.smk'
include: 'run_kmer_analysis.smk'
include: 'run_illumina_qv.smk'
include: 'run_tech_comparison.smk'

localrules: master_eval

rule master_eval:
    input:
        tech_comparison_determine_targets,
        kmer_analysis_determine_targets,
        illumina_qv_determine_targets