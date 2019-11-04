
include: 'canonical_dga.smk'
include: 'strandseq_dga_joint.smk'
include: 'strandseq_dga_split.smk'
include: 'arrow_polishing.smk'
include: 'racon_polishing.smk'
include: 'eval_known_reference.smk'
include: 'eval_variant_calls.smk'

localrules: master_results_parents, \
            create_parents_ccs_sqa_assemblies, \
            create_parents_ccs_clust_assemblies


CCS_ASSM_732_WTDBG = 'HG00732_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version'])
CCS_ASSM_732_FLYE = 'HG00732_hgsvc_pbsq2-ccs_scV{}-flye'.format(config['git_commit_version'])

CCS_ASSM_731_WTDBG = 'HG00731_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version'])
CCS_ASSM_731_FLYE = 'HG00731_hgsvc_pbsq2-ccs_scV{}-flye'.format(config['git_commit_version'])

CCS_ASSM_PARENTS = [CCS_ASSM_731_WTDBG, CCS_ASSM_731_FLYE, CCS_ASSM_732_WTDBG, CCS_ASSM_732_FLYE]


rule create_parents_ccs_sqa_assemblies:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/squashed/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                reference=[r.replace('scV{}'.format(config['git_commit_version']), 'sqa') for r in CCS_ASSM_PARENTS]
                )

rule create_parents_ccs_clust_assemblies:
    input:
        expand('output/evaluation/quastlg_busco/{known_ref}-{genemodel}/reference_assembly/clustered/{sts_reads}/{reference}/report.pdf',
                known_ref=['GRCh38_GCA_p13'],
                genemodel=['GRCh38_GENCODEv31_basic'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                reference=CCS_ASSM_PARENTS
                )


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
        rules.create_parents_ccs_sqa_assemblies.input,
        rules.create_parents_ccs_clust_assemblies.input


#         expand('input/fastq/strand-seq/{individual}_{bioproject}/requests',
#                 individual=['HG00732', 'HG00731'],
#                 bioproject=['PRJEB12849']),
#
#         expand('output/plotting/statistics/fastq_input/{filename}.stats.pdf',
#                 filename=['HG00731_hgsvc_pbsq2-ccs_1000', 'HG00732_hgsvc_pbsq2-ccs_1000']),
#
#         expand('output/evaluation/known_reference/quastlg_busco/{reference}.{known_ref}.{genemodel}/report.pdf',
#                 reference=[
#                     'HG00731_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version']),
#                     'HG00732_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version']),
#                     'HG00731_hgsvc_pbsq2-ccs_sqa-wtdbg',
#                     'HG00732_hgsvc_pbsq2-ccs_sqa-wtdbg',
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
#                 reference=['HG00731_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version'])],
#                 vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
#                 hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
#                 sts_reads=['HG00731_1kg_il25k-npe_sseq'],
#                 pol_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
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
#                 reference=['HG00732_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version'])],
#                 vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
#                 hap_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
#                 sts_reads=['HG00732_1kg_il25k-npe_sseq'],
#                 pol_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
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
#                 reference=['HG00731_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version'])],
#                 vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
#                 hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
#                 sts_reads=['HG00731_1kg_il25k-npe_sseq'],
#                 hap=['h1-un', 'h2-un', 'h1', 'h2']
#                 ),
#
#         expand('output/plotting/statistics/fastq_haplosplit/{approach}_joint/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.stats.pdf',
#                 approach=['strandseq'],
#                 var_caller=['freebayes'],
#                 gq=config['filter_vcf_gq']['freebayes'],
#                 dp=config['filter_vcf_dp']['freebayes'],
#                 reference=['HG00732_hgsvc_pbsq2-ccs_scV{}-wtdbg'.format(config['git_commit_version'])],
#                 vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
#                 hap_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
#                 sts_reads=['HG00732_1kg_il25k-npe_sseq'],
#                 hap=['h1-un', 'h2-un', 'h1', 'h2']
#                 ),
