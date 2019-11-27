
include: 'aux_utilities.smk'
include: 'canonical_dga.smk'
include: 'strandseq_dga_joint.smk'
include: 'strandseq_dga_split.smk'
include: 'integrative_phasing.smk'
include: 'arrow_polishing.smk'
include: 'racon_polishing.smk'
include: 'eval_known_reference.smk'
include: 'eval_variant_calls.smk'

localrules: master_results_child, \
            create_child_clr_nhr_assemblies, \
            create_child_ccs_nhr_assemblies, \
            create_child_clr_clust_assemblies, \
            create_child_ccs_clust_assemblies, \
            child_733_pure_ccs_variant_calling, \
            child_733_pure_clr_variant_calling, \
            child_733_pure_ccs_split_sdga, \
            child_733_pure_clr_split_sdga, \
            child_733_mixed_clr_split_sdga


USE_SINGULARITY = bool(config['use_singularity'])


CLR_ASSM_733_WTDBG = 'HG00733_sra_pbsq1-clr_1000_scV{}-wtdbg'.format(config['git_commit_version'])
CLR_ASSM_733_FLYE = 'HG00733_sra_pbsq1-clr_1000_scV{}-flye'.format(config['git_commit_version'])
CLR_ASSM_733 = [CLR_ASSM_733_WTDBG, CLR_ASSM_733_FLYE]


CCS_ASSM_733_WTDBG = 'HG00733_hgsvc_pbsq2-ccs_1000_scV{}-wtdbg'.format(config['git_commit_version'])
CCS_ASSM_733_FLYE = 'HG00733_hgsvc_pbsq2-ccs_1000_scV{}-flye'.format(config['git_commit_version'])
CCS_ASSM_733_PEREG = 'HG00733_hgsvc_pbsq2-ccs_1000_scV{}-pereg'.format(config['git_commit_version'])

if USE_SINGULARITY:
    CCS_ASSM_733 = [CCS_ASSM_733_PEREG]
    CCS_VAR_CALLER = 'deepvar'
else:
    CCS_ASSM_733 = [CCS_ASSM_733_WTDBG, CCS_ASSM_733_FLYE]
    CCS_VAR_CALLER = 'freebayes'


if not USE_SINGULARITY:
    rule create_child_clr_nhr_assemblies:
        input:
            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/non-hap-res/{reference}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    reference=[r.replace('scV{}'.format(config['git_commit_version']), 'nhr') for r in CLR_ASSM_733]
                    )


    rule create_child_clr_clust_assemblies:
        input:
            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/clustered/{sts_reads}/{reference}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    reference=CLR_ASSM_733
                    )

    rule child_733_pure_clr_variant_calling:
        input:
            expand('output/statistics/variant_calls/{var_caller}/{reference}/{sts_reads}/{vc_reads}.snv.QUAL{qual}.GQ{gq}.vcf.stats',
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733,
                    vc_reads=['HG00733_sra_pbsq1-clr_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                   ),

            expand('output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.{ext}',
                    var_caller=['longshot'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733,
                    vc_reads=['HG00733_sra_pbsq1-clr_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                    ext=['vcf.stats', 'stats.tsv', 'stats.txt']
                   ),

else:
    rule create_child_clr_nhr_assemblies:
        input:
            rules.no_singularity_mock_output.output

    rule create_child_clr_clust_assemblies:
        input:
            rules.no_singularity_mock_output.output

    rule child_733_pure_clr_variant_calling:
        input:
            rules.no_singularity_mock_output.output


rule create_child_ccs_nhr_assemblies:
        input:
            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/non-hap-res/{reference}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    reference=[r.replace('scV{}'.format(config['git_commit_version']), 'nhr') for r in CCS_ASSM_733]
                    )


rule create_child_ccs_clust_assemblies:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/clustered/{sts_reads}/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                reference=CCS_ASSM_733
                )


rule child_733_pure_ccs_variant_calling:
    input:
        expand('output/statistics/variant_calls/{var_caller}/{reference}/{sts_reads}/{vc_reads}.snv.QUAL{qual}.GQ{gq}.vcf.stats',
                var_caller=[CCS_VAR_CALLER],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_733,
                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
               ),

        expand('output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.{ext}',
                var_caller=[CCS_VAR_CALLER],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_733,
                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                ext=['vcf.stats', 'stats.tsv', 'stats.txt']
               ),


if not USE_SINGULARITY:
    rule child_733_pure_clr_split_sdga:
        input:
            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
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

            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
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

            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{pol_pass}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
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

            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{pol_pass}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
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

        priority: 250


    rule child_733_pure_ccs_split_sdga:
        input:
            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    var_caller=['freebayes'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CCS_ASSM_733_WTDBG,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    assembler=['wtdbg'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    var_caller=['freebayes'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CCS_ASSM_733_FLYE,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    assembler=['flye'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),


            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{pol_pass}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    var_caller=['freebayes'],
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

            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{pol_pass}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    var_caller=['freebayes'],
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



    rule child_733_mixed_clr_split_sdga:
        input:
            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    var_caller=['freebayes'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733_FLYE,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                    assembler=['flye'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    var_caller=['freebayes'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=CLR_ASSM_733_WTDBG,
                    vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                    assembler=['wtdbg'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{pol_pass}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    var_caller=['freebayes'],
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

            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{pol_pass}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    var_caller=['freebayes'],
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
    rule child_733_pure_clr_split_sdga:
        input:
            rules.no_singularity_mock_output.output

    rule child_733_pure_ccs_split_sdga:
        input:
            rules.no_singularity_mock_output.output

    rule child_733_mixed_clr_split_sdga:
        input:
            rules.no_singularity_mock_output.output



rule master_results_child:
    """
    Available datasets for child individuals:
    HG00733:
        - HG00733_hgsvc_pbsq2-ccs_1000.fastq.gz
        - HG00733_sra_pbsq1-clr_1000.fastq.gz
        - HG00733_hpg_ontpm-ul_1000.fastq.gz
        - HG00733_1kg_il25k-npe_sseq.fastq.gz
    HG002:
        - HG002_pbio_pbsq2-ccs_1000.fastq.gz
        - HG002_v19_pbsq2-ccs_1925.fastq.gz
        - HG002_v20_pbsq2-ccs_1520.fastq.gz
    """
    input:
        rules.create_child_clr_nhr_assemblies.input,
        rules.create_child_ccs_nhr_assemblies.input,
        rules.create_child_clr_clust_assemblies.input,
        rules.create_child_ccs_clust_assemblies.input,
        rules.child_733_pure_clr_variant_calling.input,
        rules.child_733_pure_ccs_variant_calling.input,
        rules.child_733_pure_ccs_split_sdga.input,
        rules.child_733_pure_clr_split_sdga.input,
        rules.child_733_mixed_clr_split_sdga.input
