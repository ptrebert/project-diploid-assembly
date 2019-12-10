
include: '../result_vars.smk'

localrules: run_chs_child,
            chs_child_clr_nhr_assemblies,
            chs_child_ccs_nhr_assemblies,
            chs_child_clr_clust_assemblies,
            chs_child_ccs_clust_assemblies,
            chs_child_clr_variant_calling,
            chs_child_ccs_variant_calling,
            chs_child_ccs_split_sdga,
            chs_child_clr_split_sdga_wtdbg,
            chs_child_clr_split_sdga_flye,
            chs_child_clrccs_split_sdga

USE_SINGULARITY = bool(config['use_singularity'])

CLR_ASSM_514_WTDBG = 'HG00514_hgsvc_pbsq2-clr_1000_scV{}-wtdbg'.format(config['git_commit_version'])
CLR_ASSM_514_FLYE = 'HG00514_hgsvc_pbsq2-clr_1000_scV{}-flye'.format(config['git_commit_version'])
CLR_ASSM_514_SHASTA = 'HG00514_hgsvc_pbsq2-clr_1000_scV{}-shasta'.format(config['git_commit_version'])
CLR_ASSM_514 = [CLR_ASSM_514_WTDBG, CLR_ASSM_514_FLYE, CLR_ASSM_514_SHASTA]

CCS_ASSM_514_WTDBG = 'HG00514_hgsvc_pbsq2-ccs_1000_scV{}-wtdbg'.format(config['git_commit_version'])
CCS_ASSM_514_FLYE = 'HG00514_hgsvc_pbsq2-ccs_1000_scV{}-flye'.format(config['git_commit_version'])
CCS_ASSM_514_PEREG = 'HG00514_hgsvc_pbsq2-ccs_1000_scV{}-pereg'.format(config['git_commit_version'])

if USE_SINGULARITY:
    CCS_ASSM_514 = [CCS_ASSM_514_PEREG]
else:
    CCS_ASSM_514 = [CCS_ASSM_514_WTDBG, CCS_ASSM_514_FLYE]

if not USE_SINGULARITY:
    rule chs_child_clr_nhr_assemblies:
        input:
            expand(RESULT_NHR_ASSEMBLY_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    reference=[r.replace('scV{}'.format(config['git_commit_version']), 'nhr') for r in CLR_ASSM_514]
                    )


    rule chs_child_clr_clust_assemblies:
        input:
            expand(RESULT_CLUSTERED_ASSEMBLY_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    reference=CLR_ASSM_514
                    )

    rule chs_child_clr_variant_calling:
        input:
            expand(RESULT_VARIANT_CALL_STATS,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_514,
                    vc_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                   ),

            expand(RESULT_PHASING_STATS,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_514,
                    vc_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-clr_1000']
                   ),

else:
    rule chs_child_clr_nhr_assemblies:
        input:
            rules.no_singularity_mock_output.output

    rule chs_child_clr_clust_assemblies:
        input:
            rules.no_singularity_mock_output.output

    rule chs_child_clr_variant_calling:
        input:
            rules.no_singularity_mock_output.output


rule chs_child_ccs_nhr_assemblies:
        input:
            expand(RESULT_NHR_ASSEMBLY_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    reference=[r.replace('scV{}'.format(config['git_commit_version']), 'nhr') for r in CCS_ASSM_514]
                    )


rule chs_child_ccs_clust_assemblies:
    input:
        expand(RESULT_CLUSTERED_ASSEMBLY_QUAST_REPORT,
                known_ref=KNOWN_REF,
                genemodel=GENEMODEL,
                sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                reference=CCS_ASSM_514
                )


rule chs_child_ccs_variant_calling:
    input:
        expand(RESULT_VARIANT_CALL_STATS,
                var_caller=CCS_VAR_CALLER,
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_514,
                vc_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00514_1kg_il25k-npe_sseq'],
               ),

        expand(RESULT_PHASING_STATS,
                var_caller=CCS_VAR_CALLER,
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_514,
                vc_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                hap_reads=['HG00514_hgsvc_pbsq2-ccs_1000']
               ),


if not USE_SINGULARITY:
    rule chs_child_clr_split_sdga_wtdbg:
        input:
            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_514_WTDBG,
                    vc_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    assembler=['wtdbg'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_514_WTDBG,
                    vc_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    pol_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    assembler=['wtdbg'],
                    pol_pass=['arrow-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

        priority: 500


    rule chs_child_clr_split_sdga_flye:
        input:
            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_514_FLYE,
                    vc_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    assembler=['flye'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_514_FLYE,
                    vc_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    pol_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    assembler=['flye'],
                    pol_pass=['arrow-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

        priority: 500


    rule chs_child_clr_split_sdga_shasta:
        input:
            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_514_SHASTA,
                    vc_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    assembler=['shasta'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_514_SHASTA,
                    vc_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    pol_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    assembler=['shasta'],
                    pol_pass=['arrow-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

        priority: 500


    rule chs_child_ccs_split_sdga:
        input:
            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CCS_ASSM_514_WTDBG,
                    vc_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    assembler=['wtdbg'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CCS_ASSM_514_FLYE,
                    vc_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    assembler=['flye'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),


            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CCS_ASSM_514_WTDBG,
                    vc_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    pol_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
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
                    reference=CCS_ASSM_514_FLYE,
                    vc_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    pol_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    assembler=['flye'],
                    pol_pass=['racon-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

        priority: 1000


    rule chs_child_clrccs_split_sdga:
        input:
            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_514_FLYE,
                    vc_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    assembler=['flye'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_DRAFT_QUAST_REPORT ,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_514_WTDBG,
                    vc_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    assembler=['wtdbg'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand(RESULT_SPLIT_DGA_POLISHED_QUAST_REPORT,
                    known_ref=KNOWN_REF,
                    genemodel=GENEMODEL,
                    var_caller=CCS_VAR_CALLER,
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_514_WTDBG,
                    vc_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    pol_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
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
                    reference=CLR_ASSM_514_FLYE,
                    vc_reads=['HG00514_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00514_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    pol_reads=['HG00514_hgsvc_pbsq2-clr_1000'],
                    assembler=['flye'],
                    pol_pass=['arrow-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

        priority: 100

else:
    rule chs_child_clr_split_sdga_wtdbg:
        input:
            rules.no_singularity_mock_output.output

    rule chs_child_clr_split_sdga_flye:
        input:
            rules.no_singularity_mock_output.output

    rule chs_child_ccs_split_sdga:
        input:
            rules.no_singularity_mock_output.output

    rule chs_child_clrccs_split_sdga:
        input:
            rules.no_singularity_mock_output.output


rule run_chs_child:
    input:
        rules.chs_child_clr_nhr_assemblies.input,
        rules.chs_child_ccs_nhr_assemblies.input,
        rules.chs_child_clr_clust_assemblies.input,
        rules.chs_child_ccs_clust_assemblies.input,
        rules.chs_child_clr_variant_calling.input,
        rules.chs_child_ccs_variant_calling.input,
        rules.chs_child_ccs_split_sdga.input,
        rules.chs_child_clr_split_sdga_wtdbg.input,
        rules.chs_child_clr_split_sdga_flye.input,
        rules.chs_child_clrccs_split_sdga.input