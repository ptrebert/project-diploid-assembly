
include: 'canonical_dga.smk'
include: 'strandseq_dga_joint.smk'
include: 'strandseq_dga_split.smk'
include: 'arrow_polishing.smk'
include: 'racon_polishing.smk'
include: 'collect_statistics.smk'
include: 'eval_known_reference.smk'
include: 'eval_variant_calls.smk'

localrules: master_results_parents, \
            parents_ccs_sqa_assemblies, \
            parents_ccs_clust_assemblies, \
            parent_731_pure_ccs_split_sdga, \
            parent_732_pure_ccs_split_sdga


CCS_ASSM_732_WTDBG = 'HG00732_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version'])
CCS_ASSM_732_FLYE = 'HG00732_hgsvc_pbsq2-ccs_scV{}-flye'.format(config['git_commit_version'])

CCS_ASSM_731_WTDBG = 'HG00731_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version'])
CCS_ASSM_731_FLYE = 'HG00731_hgsvc_pbsq2-ccs_scV{}-flye'.format(config['git_commit_version'])

CCS_ASSM_PARENTS = [CCS_ASSM_731_WTDBG, CCS_ASSM_731_FLYE, CCS_ASSM_732_WTDBG, CCS_ASSM_732_FLYE]


rule parents_ccs_sqa_assemblies:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/squashed/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                reference=[r.replace('scV{}'.format(config['git_commit_version']), 'sqa') for r in CCS_ASSM_PARENTS]
                )

rule parents_ccs_clust_assemblies:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/clustered/{sts_reads}/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                sts_reads=['HG00732_1kg_il25k-npe_sseq'],
                reference=[CCS_ASSM_732_WTDBG, CCS_ASSM_732_FLYE]
                ),

        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/clustered/{sts_reads}/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                reference=[CCS_ASSM_731_WTDBG, CCS_ASSM_731_FLYE]
                )


rule parent_731_pure_ccs_split_sdga:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                var_caller=['freebayes'],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=[CCS_ASSM_731_FLYE],
                vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                assembler=['flye'],
                hap=['h1-un', 'h2-un', 'h1', 'h2']
               ),

        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                var_caller=['freebayes'],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=[CCS_ASSM_731_WTDBG],
                vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                assembler=['wtdbg'],
                hap=['h1-un', 'h2-un', 'h1', 'h2']
               ),

        expand('output/statistics/variant_calls/{var_caller}/{reference}/{sts_reads}/{vc_reads}.snps.QUAL{qual}.GQ{gq}.vcf.stats',
                var_caller=['freebayes'],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=[CCS_ASSM_731_FLYE, CCS_ASSM_731_WTDBG],
                vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
               ),

        expand('output/statistics/variant_calls/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.snps.phased.vcf.stats',
                var_caller=['freebayes'],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=[CCS_ASSM_731_FLYE, CCS_ASSM_731_WTDBG],
                vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000']
               ),

        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{pol_pass}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                var_caller=['freebayes'],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_731_WTDBG,
                vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                pol_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
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
                reference=CCS_ASSM_731_FLYE,
                vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                pol_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                assembler=['flye'],
                pol_pass=['racon-p1'],
                hap=['h1-un', 'h2-un', 'h1', 'h2']
               ),

    priority: 1000


rule parent_732_pure_ccs_split_sdga:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                var_caller=['freebayes'],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=[CCS_ASSM_732_FLYE],
                vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00732_1kg_il25k-npe_sseq'],
                hap_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                assembler=['flye'],
                hap=['h1-un', 'h2-un', 'h1', 'h2']
               ),

        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                var_caller=['freebayes'],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=[CCS_ASSM_732_WTDBG],
                vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00732_1kg_il25k-npe_sseq'],
                hap_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                assembler=['wtdbg'],
                hap=['h1-un', 'h2-un', 'h1', 'h2']
               ),

        expand('output/statistics/variant_calls/{var_caller}/{reference}/{sts_reads}/{vc_reads}.snps.QUAL{qual}.GQ{gq}.vcf.stats',
                var_caller=['freebayes'],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=[CCS_ASSM_732_FLYE, CCS_ASSM_732_WTDBG],
                vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00732_1kg_il25k-npe_sseq'],
               ),

        expand('output/statistics/variant_calls/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.snps.phased.vcf.stats',
                var_caller=['freebayes'],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=[CCS_ASSM_732_FLYE, CCS_ASSM_732_WTDBG],
                vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00732_1kg_il25k-npe_sseq'],
                hap_reads=['HG00732_hgsvc_pbsq2-ccs_1000']
               ),

        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{pol_pass}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                var_caller=['freebayes'],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_732_WTDBG,
                vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00732_1kg_il25k-npe_sseq'],
                hap_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                pol_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
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
                reference=CCS_ASSM_732_FLYE,
                vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00732_1kg_il25k-npe_sseq'],
                hap_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                pol_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                assembler=['flye'],
                pol_pass=['racon-p1'],
                hap=['h1-un', 'h2-un', 'h1', 'h2']
               ),

    priority: 1000


rule master_results_parents:
    """
    Available datasets for parent individuals:
    HG00731:
        - HG00731_hgsvc_pbsq2-ccs_1000.fastq.gz
        - HG00731_1kg_il25k-npe_sseq.fastq.gz
    HG00732:
        - HG00732_1kg_il25k-npe_sseq.fastq.gz
        - HG00732_hgsvc_pbsq2-ccs_1000.fastq.gz
    """
    input:
        rules.parents_ccs_sqa_assemblies.input,
        rules.parents_ccs_clust_assemblies.input,
        rules.parent_731_pure_ccs_split_sdga.input,
        rules.parent_732_pure_ccs_split_sdga.input
