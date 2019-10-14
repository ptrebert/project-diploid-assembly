
include: 'eval_known_reference.smk'
include: 'eval_variant_calls.smk'

localrules: master_results_parents

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
#        # DONE
#        expand('output/evaluation/known_reference/quastlg_busco/{approach}.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                approach=['strandseq'],
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00731_hgsvc_pbsq2-ccs_clustV2-100kb'],
#                vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
#                sts_reads=['HG00731_1kg_il25k-npe_sseq'],
#                hap=['h1', 'h2', 'h1-un', 'h2-un'],
#                known_ref=['GRCh38_GCA_p13'],
#                genemodel=['GRCh38_GENCODEv31_basic']
#                ),
#
#        # DONE
#        expand('output/evaluation/known_reference/quastlg_busco/{approach}.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                approach=['strandseq'],
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00732_hgsvc_pbsq2-ccs_clustV2-100kb'],
#                vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
#                sts_reads=['HG00732_1kg_il25k-npe_sseq'],
#                hap=['h1', 'h2', 'h1-un', 'h2-un'],
#                known_ref=['GRCh38_GCA_p13'],
#                genemodel=['GRCh38_GENCODEv31_basic']
#                ),
