
KNOWN_REF = ['GRCh38_GCA_p13']
GENEMODEL = ['GRCh38_GENCODEv31_basic']

USE_SINGULARITY = bool(config['use_singularity'])

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