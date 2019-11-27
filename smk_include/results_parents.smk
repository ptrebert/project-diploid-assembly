
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


USE_SINGULARITY = bool(config['use_singularity'])


CCS_ASSM_732_WTDBG = 'HG00732_hgsvc_pbsq2-ccs_1000_scV{}-wtdbg'.format(config['git_commit_version'])
CCS_ASSM_732_FLYE = 'HG00732_hgsvc_pbsq2-ccs_1000_scV{}-flye'.format(config['git_commit_version'])
CCS_ASSM_732_PEREG = 'HG00732_hgsvc_pbsq2-ccs_1000_scV{}-pereg'.format(config['git_commit_version'])
CCS_ASSM_732 = [ CCS_ASSM_732_WTDBG , CCS_ASSM_732_FLYE ]

if USE_SINGULARITY:
    CCS_ASSM_732 = [CCS_ASSM_732_PEREG]


CCS_ASSM_731_WTDBG = 'HG00731_hgsvc_pbsq2-ccs_1000_scV{}-wtdbg'.format(config['git_commit_version'])
CCS_ASSM_731_FLYE = 'HG00731_hgsvc_pbsq2-ccs_1000_scV{}-flye'.format(config['git_commit_version'])
CCS_ASSM_731_PEREG = 'HG00731_hgsvc_pbsq2-ccs_1000_scV{}-pereg'.format(config['git_commit_version'])
CCS_ASSM_731 = [ CCS_ASSM_731_WTDBG , CCS_ASSM_731_FLYE ]

if USE_SINGULARITY:
    CCS_ASSM_731 = [ CCS_ASSM_731_PEREG ]


CCS_ASSM_PARENTS = [CCS_ASSM_731_WTDBG, CCS_ASSM_731_FLYE, CCS_ASSM_732_WTDBG, CCS_ASSM_732_FLYE]
CCS_VAR_CALLER = 'freebayes'

if USE_SINGULARITY:
    CCS_ASSM_PARENTS = [ CCS_ASSM_731_PEREG , CCS_ASSM_732_PEREG ]
    CCS_VAR_CALLER = 'deepvar'


rule parents_ccs_nhr_assemblies:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/non-hap-res/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                reference=[r.replace('scV{}'.format(config['git_commit_version']), 'nhr') for r in CCS_ASSM_PARENTS]
                )

rule parents_ccs_clust_assemblies:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/clustered/{sts_reads}/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                sts_reads=['HG00732_1kg_il25k-npe_sseq'],
                reference=CCS_ASSM_PARENTS
                ),

        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/clustered/{sts_reads}/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                reference=CCS_ASSM_PARENTS
                )


rule parent_731_pure_ccs_variant_calling:
    input:
        expand('output/statistics/variant_calls/{var_caller}/{reference}/{sts_reads}/{vc_reads}.snv.QUAL{qual}.GQ{gq}.vcf.stats',
            var_caller=[CCS_VAR_CALLER],
            qual=config['filter_vcf_qual'],
            gq=config['filter_vcf_gq'],
            reference=CCS_ASSM_731,
            vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
            sts_reads=['HG00731_1kg_il25k-npe_sseq'],
           ),

        expand('output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.{ext}',
                var_caller=[CCS_VAR_CALLER],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_731,
                vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
                hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
                ext=['vcf.stats', 'stats.tsv', 'stats.txt']
               ),


rule parent_732_pure_ccs_variant_calling:
    input:
        expand('output/statistics/variant_calls/{var_caller}/{reference}/{sts_reads}/{vc_reads}.snv.QUAL{qual}.GQ{gq}.vcf.stats',
            var_caller=[CCS_VAR_CALLER],
            qual=config['filter_vcf_qual'],
            gq=config['filter_vcf_gq'],
            reference=CCS_ASSM_732,
            vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
            sts_reads=['HG00732_1kg_il25k-npe_sseq'],
           ),

        expand('output/statistics/phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.{ext}',
                var_caller=[CCS_VAR_CALLER],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_732,
                vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00732_1kg_il25k-npe_sseq'],
                hap_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                ext=['vcf.stats', 'stats.tsv', 'stats.txt']
               ),


if not USE_SINGULARITY:
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

            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{pol_pass}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    var_caller=['freebayes'],
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

            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{pol_pass}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    var_caller=['freebayes'],
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

            expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{pol_pass}/report.pdf',
                    known_ref=['GRCh38_GCA_p13'],
                    genemodel=['GRCh38_GENCODEv31_basic'],
                    var_caller=['freebayes'],
                    qual=config['filter_vcf_qual'],
                    gq=config['filter_vcf_gq'],
                    reference=[CCS_ASSM_732_WTDBG],
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
                    reference=[CCS_ASSM_732_FLYE],
                    vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                    sts_reads=['HG00732_1kg_il25k-npe_sseq'],
                    hap_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                    pol_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
                    assembler=['flye'],
                    pol_pass=['racon-p1'],
                    hap=['h1-un', 'h2-un', 'h1', 'h2']
                   ),

        priority: 1000

else:
    rule parent_731_pure_ccs_split_sdga:
        input:
            rules.no_singularity_mock_output.output

    rule parent_732_pure_ccs_split_sdga:
        input:
            rules.no_singularity_mock_output.output


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
        rules.parents_ccs_nhr_assemblies.input,
        rules.parents_ccs_clust_assemblies.input,
        rules.parent_731_pure_ccs_variant_calling.input,
        rules.parent_732_pure_ccs_variant_calling.input,
        rules.parent_731_pure_ccs_split_sdga.input,
        rules.parent_732_pure_ccs_split_sdga.input
