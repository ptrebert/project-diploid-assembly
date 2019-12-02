
include: '../result_vars.smk'

localrules: run_pur_child,
            pur_father_ccs_nhr_assemblies,
            pur_father_ccs_clust_assemblies,
            pur_father_ccs_variant_calling,
            pur_father_ccs_split_sdga,

USE_SINGULARITY = bool(config['use_singularity'])

CCS_ASSM_731_WTDBG = 'HG00731_hgsvc_pbsq2-ccs_1000_scV{}-wtdbg'.format(config['git_commit_version'])
CCS_ASSM_731_FLYE = 'HG00731_hgsvc_pbsq2-ccs_1000_scV{}-flye'.format(config['git_commit_version'])
CCS_ASSM_731_PEREG = 'HG00731_hgsvc_pbsq2-ccs_1000_scV{}-pereg'.format(config['git_commit_version'])
CCS_ASSM_731 = [ CCS_ASSM_731_WTDBG , CCS_ASSM_731_FLYE ]

if USE_SINGULARITY:
    CCS_ASSM_731 = [ CCS_ASSM_731_PEREG ]

rule pur_father_ccs_nhr_assemblies:
    input:
        expand(RESULT_NHR_ASSEMBLY_QUAST_REPORT,
                known_ref=KNOWN_REF,
                genemodel=GENEMODEL,
                reference=[r.replace('scV{}'.format(config['git_commit_version']), 'nhr') for r in CCS_ASSM_731]
                )

rule pur_father_ccs_clust_assemblies:
    input:
        expand(RESULT_CLUSTERED_ASSEMBLY_QUAST_REPORT,
                known_ref=KNOWN_REF,
                genemodel=GENEMODEL,
                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                reference=CCS_ASSM_731
                )

rule pur_father_ccs_variant_calling:
    input:
        expand(RESULT_VARIANT_CALL_STATS,
                var_caller=CCS_VAR_CALLER,
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_731,
                vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
               ),

        expand(RESULT_PHASING_STATS,
                var_caller=CCS_VAR_CALLER,
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_731,
                vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
               ),


if not USE_SINGULARITY:
    rule pur_father_ccs_split_sdga:
        input:
            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=[CCS_ASSM_731_FLYE],
                    vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                    assembler=['flye'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=[CCS_ASSM_731_WTDBG],
                    vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                    assembler=['wtdbg'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=[CCS_ASSM_731_WTDBG],
                    vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                    pol_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                    assembler=['wtdbg'],
                    pol_pass=['racon-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=[CCS_ASSM_731_FLYE],
                    vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                    pol_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                    assembler=['flye'],
                    pol_pass=['racon-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

        priority: 1000


else:
    rule pur_father_ccs_split_sdga:
        input:
            rules.no_singularity_mock_output.output


rule run_pur_father:
    input:
        rules.pur_father_ccs_nhr_assemblies.input,
        rules.pur_father_ccs_clust_assemblies.input,
        rules.pur_father_ccs_variant_calling.input,
        rules.pur_father_ccs_split_sdga.input,