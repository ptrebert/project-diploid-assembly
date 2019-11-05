
include: 'canonical_dga.smk'
include: 'strandseq_dga_joint.smk'
include: 'strandseq_dga_split.smk'
include: 'integrative_phasing.smk'
include: 'arrow_polishing.smk'
include: 'racon_polishing.smk'
include: 'eval_known_reference.smk'
include: 'eval_variant_calls.smk'

localrules: master_results_child, \
            create_child_clr_sqa_assemblies, \
            create_child_ccs_sqa_assemblies, \
            create_child_clr_clust_assemblies, \
            create_child_ccs_clust_assemblies


CLR_ASSM_733_WTDBG = 'HG00733_sra_pbsq1-clr_scV{}-wtdbg'.format(config['git_commit_version'])
CLR_ASSM_733_FLYE = 'HG00733_sra_pbsq1-clr_scV{}-flye'.format(config['git_commit_version'])
CLR_ASSM_733 = [CLR_ASSM_733_WTDBG, CLR_ASSM_733_FLYE]


CCS_ASSM_733_WTDBG = 'HG00733_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version'])
CCS_ASSM_733_FLYE = 'HG00733_hgsvc_pbsq2-ccs_scV{}-flye'.format(config['git_commit_version'])
CCS_ASSM_733 = [CCS_ASSM_733_WTDBG, CCS_ASSM_733_FLYE]


rule create_child_clr_sqa_assemblies:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/squashed/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                reference=[r.replace('scV{}'.format(config['git_commit_version']), 'sqa') for r in CLR_ASSM_733]
                )

rule create_child_ccs_sqa_assemblies:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/squashed/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                reference=[r.replace('scV{}'.format(config['git_commit_version']), 'sqa') for r in CCS_ASSM_733]
                )

rule create_child_clr_clust_assemblies:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/clustered/{sts_reads}/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                reference=CLR_ASSM_733
                )

rule create_child_ccs_clust_assemblies:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/clustered/{sts_reads}/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                reference=CCS_ASSM_733
                )


rule perform_child_pure_clr_sdga:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                var_caller=['longshot'],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CLR_ASSM_733,
                vc_reads=['HG00733_sra_pbsq1-clr_1000'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                assembler=['flye', 'wtdbg'],
                hap=['h1-un', 'h2-un', 'h1', 'h2']
               )


rule perform_child_pure_ccs_sdga:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                var_caller=['freebayes'],
                qual=config['filter_vcf_qual'],
                gq=config['filter_vcf_gq'],
                reference=CCS_ASSM_733,
                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
                assembler=['flye', 'wtdbg'],
                hap=['h1-un', 'h2-un', 'h1', 'h2']
               )


rule perform_child_mixed_clr_sdga:
    input:
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
                assembler=['wtdbg', 'flye'],
                hap=['h1-un', 'h2-un', 'h1', 'h2']
               )


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
        rules.create_child_clr_sqa_assemblies.input,
        rules.create_child_ccs_sqa_assemblies.input,
        rules.create_child_clr_clust_assemblies.input,
        rules.create_child_ccs_clust_assemblies.input,
        rules.perform_child_pure_clr_sdga.input,
        rules.perform_child_pure_ccs_sdga.input,
        rules.perform_child_mixed_clr_sdga.input






#         expand('input/fastq/strand-seq/{individual}_{bioproject}/requests',
#                 individual=['HG00733'],
#                 bioproject=['PRJEB12849']),
#
#         expand('output/plotting/statistics/fastq_input/{filename}.stats.pdf',
#                 filename=['HG00733_sra_pbsq1-clr_1000', 'HG00733_hgsvc_pbsq2-ccs_1000']),
#
#         expand('output/evaluation/known_reference/quastlg_busco/{reference}.{known_ref}.{genemodel}/report.pdf',
#                 reference=[
#                     'HG00733_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version']),
#                     'HG00733_sra_pbsq1-clr_scV{}-wtdbg'.format(config['git_commit_version']),
#                     'HG00733_hgsvc_pbsq2-ccs_sqa-wtdbg',
#                     'HG00733_sra_pbsq1-clr_sqa-wtdbg',
#                 ],
#                 known_ref=['GRCh38_GCA_p13'],
#                 genemodel=['GRCh38_GENCODEv31_basic']
#                 ),
#
#         expand('output/evaluation/known_reference/quastlg_busco/strandseq.split.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{assembler}.{pol_reads}.{polisher}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                 var_caller=['freebayes'],
#                 assembler=['wtdbg'],
#                 gq=config['filter_vcf_gq']['freebayes'],
#                 dp=config['filter_vcf_dp']['freebayes'],
#                 reference=['HG00733_sra_pbsq1-clr_scV{}-wtdbg'.format(config['git_commit_version'])],
#                 vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                 hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                 sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                 pol_reads=['HG00733_sra_pbsq1-clr_1000'],
#                 polisher=['arrow-p1'],
#                 hap=['h1-un', 'h2-un', 'h1', 'h2'],
#                 known_ref=['GRCh38_GCA_p13'],
#                 genemodel=['GRCh38_GENCODEv31_basic']
#                 ),
#
#         expand('output/evaluation/known_reference/quastlg_busco/strandseq.split.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{assembler}.{pol_reads}.{polisher}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                 var_caller=['freebayes'],
#                 assembler=['wtdbg'],
#                 gq=config['filter_vcf_gq']['freebayes'],
#                 dp=config['filter_vcf_dp']['freebayes'],
#                 reference=['HG00733_sra_pbsq1-clr_scV{}-wtdbg'.format(config['git_commit_version'])],
#                 vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                 hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                 sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                 pol_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                 polisher=['racon-p1'],
#                 hap=['h1-un', 'h2-un', 'h1', 'h2'],
#                 known_ref=['GRCh38_GCA_p13'],
#                 genemodel=['GRCh38_GENCODEv31_basic']
#                 ),
#
#         expand('output/evaluation/known_reference/quastlg_busco/strandseq.split.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{assembler}.{pol_reads}.{polisher}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                 var_caller=['freebayes'],
#                 assembler=['wtdbg'],
#                 gq=config['filter_vcf_gq']['freebayes'],
#                 dp=config['filter_vcf_dp']['freebayes'],
#                 reference=['HG00733_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version'])],
#                 vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                 hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                 sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                 pol_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                 polisher=['racon-p1'],
#                 hap=['h1-un', 'h2-un', 'h1', 'h2'],
#                 known_ref=['GRCh38_GCA_p13'],
#                 genemodel=['GRCh38_GENCODEv31_basic']
#                 ),
#
#         expand('output/plotting/statistics/fastq_haplosplit/{approach}_joint/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.stats.pdf',
#                 approach=['strandseq'],
#                 var_caller=['freebayes'],
#                 gq=config['filter_vcf_gq']['freebayes'],
#                 dp=config['filter_vcf_dp']['freebayes'],
#                 reference=['HG00733_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version'])],
#                 vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                 hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                 sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                 hap=['h1-un', 'h2-un', 'h1', 'h2']
#                 ),
#
#         expand('output/plotting/statistics/fastq_haplosplit/{approach}_joint/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.stats.pdf',
#                 approach=['strandseq'],
#                 var_caller=['freebayes'],
#                 gq=config['filter_vcf_gq']['freebayes'],
#                 dp=config['filter_vcf_dp']['freebayes'],
#                 reference=['HG00733_sra_pbsq1-clr_scV{}-wtdbg'.format(config['git_commit_version'])],
#                 vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                 hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                 sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                 hap=['h1-un', 'h2-un', 'h1', 'h2']
#                )
