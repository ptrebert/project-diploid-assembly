
KNOWN_REF = ['GRCh38_GCA_p13']
GENEMODEL = ['GRCh38_GENCODEv31_basic']

USE_SINGULARITY = False

if USE_SINGULARITY:
    CCS_VAR_CALLER = ['deepvar']
else:
    CCS_VAR_CALLER = ['freebayes']

RESULT_NHR_ASSEMBLY_QUAST_REPORT = 'output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/non-hap-res/{reference}/report.pdf'
RESULT_CLUSTERED_ASSEMBLY_QUAST_REPORT = 'output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/clustered/{sts_reads}/{reference}/report.pdf'

# Includes statistics for only quality filtered variants, and for final call set
RESULT_VARIANT_CALL_STATS = 'output/statistics/variant_calls/{var_caller}/{reference}/{sts_reads}/{vc_reads}.snv.QUAL{qual}.GQ{gq}.vcf.stats'

# Includes statistics for StrandPhaseR stage, and final WhatsHap stage
RESULT_PHASING_STATS = 'output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf.stats'

RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT = 'output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-{assembler}.{hap}/report.pdf'
RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT = 'output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{pol_pass}/report.pdf'

import os

TARGET_PATHS = {
    "BUILD_NHR_ASSEMBLY": os.path.join(
        "output", "reference_assembly", "non-hap-res",
        "{hap_reads}_nhr-{nhr_assembler}.fasta"
    ),

    "BUILD_CLUSTERED_ASSEMBLY": os.path.join(
        "output", "reference_assembly", "clustered",
        "{sts_reads}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}.fasta",
    ),

    "BUILD_DRAFT_HAPLOID_ASSEMBLY": os.path.join(
        "output", "diploid_assembly",
        "strandseq_{approach}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sts_reads}",
        "draft", "haploid_assembly",
        "{hap_reads}-{hap_assembler}.{hap}.fasta"
    ),

    "BUILD_POLISHED_HAPLOID_ASSEMBLY": os.path.join(
        "output", "diploid_assembly",
        "strandseq_{approach}",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sts_reads}",
        "polishing",
        "{pol_reads}",
        "haploid_assembly",
        "{hap_reads}-{hap_assembler}.{hap}.{pol_pass}.fasta"
    ),

    "STATS_VARIANT_CALLING": os.path.join(
        "output", "statistics", "variant_calls",
        "{var_caller}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{sts_reads}",
        "{vc_reads}.snv.QUAL{filter_vcf_qual}.GQ{filter_vcf_gq}.vcf.stats"
    ),

    "STATS_INTEGRATIVE_PHASING": os.path.join(
        "output", "statistics", "phasing",
        "{var_caller}_QUAL{qual}_GQ{gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sts_reads}",
        "{hap_reads}.wh-phased.stats.tsv"
    ),

    "STATS_READ_HAPLO_TAGGING": os.path.join(
        "output", "statistics", "tag_split",
        "{var_caller}_QUAL{filter_vcf_qual}_GQ{filter_vcf_gq}",
        "{hap_reads}_scV{git_commit_version}-{nhr_assembler}",
        "{vc_reads}",
        "{sts_reads}",
        "{hap_reads}.tags.{tag_type}.tsv"
    ),
}