
include: 'prep_custom_references.smk'
include: 'run_kmer_analysis.smk'
include: 'run_illumina_qv.smk'
include: 'run_tech_comparison.smk'
include: 'run_contig_remap.smk'
include: 'run_bng_hybrids.smk'


localrules: master_eval

rule master_eval:
    input:
        tech_comparison_determine_targets,
        kmer_analysis_determine_targets,
        illumina_qv_determine_targets,
        contig_remap_determine_targets,
        bng_hybrids_determine_targets