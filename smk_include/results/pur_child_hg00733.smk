
include: '../result_vars.smk'

localrules: run_pur_child,
            pur_child_clr_nhr_assemblies,
            pur_child_ccs_nhr_assemblies,
            pur_child_clr_clust_assemblies,
            pur_child_ccs_clust_assemblies,
            pur_child_clr_variant_calling,
            pur_child_ccs_variant_calling,
            pur_child_ccs_split_sdga,
            pur_child_clr_split_sdga_wtdbg,
            pur_child_clr_split_sdga_flye,
            pur_child_clrccs_split_sdga

USE_SINGULARITY = bool(config['use_singularity'])

CLR_ASSM_733_WTDBG = 'HG00733_sra_pbsq1-clr_1000_scV{}-wtdbg'.format(config['git_commit_version'])
CLR_ASSM_733_FLYE = 'HG00733_sra_pbsq1-clr_1000_scV{}-flye'.format(config['git_commit_version'])
CLR_ASSM_733 = [CLR_ASSM_733_WTDBG, CLR_ASSM_733_FLYE]

CCS_ASSM_733_WTDBG = 'HG00733_hgsvc_pbsq2-ccs_1000_scV{}-wtdbg'.format(config['git_commit_version'])
CCS_ASSM_733_FLYE = 'HG00733_hgsvc_pbsq2-ccs_1000_scV{}-flye'.format(config['git_commit_version'])
CCS_ASSM_733_PEREG = 'HG00733_hgsvc_pbsq2-ccs_1000_scV{}-pereg'.format(config['git_commit_version'])

if USE_SINGULARITY:
    CCS_ASSM_733 = [CCS_ASSM_733_PEREG]
else:
    CCS_ASSM_733 = [CCS_ASSM_733_WTDBG, CCS_ASSM_733_FLYE]

if not USE_SINGULARITY:
    rule pur_child_clr_nhr_assemblies:
        input:
            expand(RESULT_NHR_ASSEMBLY_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    reference=[r.replace('scV{}'.format(config['git_commit_version']), 'nhr') for r in CLR_ASSM_733]
                    )


    rule pur_child_clr_clust_assemblies:
        input:
            expand(RESULT_CLUSTERED_ASSEMBLY_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    reference=CLR_ASSM_733
                    )

    rule pur_child_clr_variant_calling:
        input:
            expand(RESULT_VARIANT_CALL_STATS,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733,
                    vc_reads=['HG00733_sra_pbsq1-clr_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                   ),

            expand(RESULT_PHASING_STATS,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733,
                    vc_reads=['HG00733_sra_pbsq1-clr_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_sra_pbsq1-clr_1000']
                   ),

else:
    rule pur_child_clr_nhr_assemblies:
        input:
            rules.no_singularity_mock_output.output

    rule pur_child_clr_clust_assemblies:
        input:
            rules.no_singularity_mock_output.output

    rule pur_child_clr_variant_calling:
        input:
            rules.no_singularity_mock_output.output


rule pur_child_ccs_nhr_assemblies:
        input:
            expand(RESULT_NHR_ASSEMBLY_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    reference=[r.replace('scV{}'.format(config['git_commit_version']), 'nhr') for r in CCS_ASSM_733]
                    )


rule pur_child_ccs_clust_assemblies:
    input:
        expand(RESULT_CLUSTERED_ASSEMBLY_QUAST_REPORT,
                known_ref=KNOWN_REF,
                genemodel=GENEMODEL,
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                reference=CCS_ASSM_733
                )


rule pur_child_ccs_variant_calling:
    input:
        expand(RESULT_VARIANT_CALL_STATS,
                var_caller=CCS_VAR_CALLER,
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_733,
                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
               ),

        expand(RESULT_PHASING_STATS,
                var_caller=CCS_VAR_CALLER,
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_733,
                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000']
               ),


if not USE_SINGULARITY:
    rule pur_child_clr_split_sdga_wtdbg:
        input:
            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733_WTDBG,
                    vc_reads=['HG00733_sra_pbsq1-clr_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                    assembler=['wtdbg'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733_WTDBG,
                    vc_reads=['HG00733_sra_pbsq1-clr_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                    pol_reads=['HG00733_sra_pbsq1-clr_1000'],
                    assembler=['wtdbg'],
                    pol_pass=['arrow-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

        priority: 250


    rule pur_child_clr_split_sdga_flye:
        input:
            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733_FLYE,
                    vc_reads=['HG00733_sra_pbsq1-clr_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                    assembler=['flye'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733_FLYE,
                    vc_reads=['HG00733_sra_pbsq1-clr_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                    pol_reads=['HG00733_sra_pbsq1-clr_1000'],
                    assembler=['flye'],
                    pol_pass=['arrow-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

        priority: 500


    rule pur_child_ccs_split_sdga:
        input:
            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CCS_ASSM_733_WTDBG,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    assembler=['wtdbg'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CCS_ASSM_733_FLYE,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    assembler=['flye'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),


            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CCS_ASSM_733_WTDBG,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    pol_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
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
                    reference=CCS_ASSM_733_FLYE,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    pol_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    assembler=['flye'],
                    pol_pass=['racon-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

        priority: 1000



    rule pur_child_clrccs_split_sdga:
        input:
            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733_FLYE,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                    assembler=['flye'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733_WTDBG,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                    assembler=['wtdbg'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733_WTDBG,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                    pol_reads=['HG00733_sra_pbsq1-clr_1000'],
                    assembler=['wtdbg'],
                    pol_pass=['arrow-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733_FLYE,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                    pol_reads=['HG00733_sra_pbsq1-clr_1000'],
                    assembler=['flye'],
                    pol_pass=['arrow-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

        priority: 100

else:
    rule pur_child_clr_split_sdga_wtdbg:
        input:
            rules.no_singularity_mock_output.output

    rule pur_child_clr_split_sdga_flye:
        input:
            rules.no_singularity_mock_output.output

    rule pur_child_clrccs_split_sdga:
        input:
            rules.no_singularity_mock_output.output

    rule pur_child_ccs_split_sdga:
        input:
            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CCS_ASSM_733_PEREG,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    assembler=['pereg'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CCS_ASSM_733_PEREG,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    pol_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    assembler=['pereg'],
                    pol_pass=['racon-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

        priority: 1000


rule run_pur_child:
    input:
        rules.pur_child_clr_nhr_assemblies.input,
        rules.pur_child_ccs_nhr_assemblies.input,
        rules.pur_child_clr_clust_assemblies.input,
        rules.pur_child_ccs_clust_assemblies.input,
        rules.pur_child_clr_variant_calling.input,
        rules.pur_child_ccs_variant_calling.input,
        rules.pur_child_ccs_split_sdga.input,
        rules.pur_child_clr_split_sdga_wtdbg.input,
        rules.pur_child_clr_split_sdga_flye.input,
        rules.pur_child_clrccs_split_sdga.input