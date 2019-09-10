
include: 'eval_known_reference.smk'
include: 'eval_variant_calls.smk'

localrules: master_results_child

rule master_results_child:
    input:
        expand('output/evaluation/known_reference/quastlg_busco/{custom_ref}.{known_ref}.{genemodel}/icarus.html',
                custom_ref=[
                    'HG00733_sra_pbsq1-clr_clustV1-100kb',
                    'HG00733_sra_pbsq1-clr_sqa-100kb',
                    'HG00733_hgsvc_pbsq2-ccs_sqa-100kb',
                    #'HG00733_hgsvc_pbsq2-ccs_clustV1-100kb'  # depends on SaarClust
                    ],
                known_ref=['GRCh38_1KG_decoy'],
                genemodel=['GRCh38_GENCODEv31_basic']
                ),

        expand('output/evaluation/variant_phasings/by_approach/{approach1}_vs_{approach2}/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}_{approach1}_vs_{sts_reads}.comp.tsv',
                approach1=['canonical'],
                approach2=['strandseq'],
                var_caller=['longshot'],
                gq=config['filter_vcf_gq'],
                dp=config['filter_vcf_dp'],
                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
                vc_reads=['HG00733_sra_pbsq1-clr_1000'],
                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq']
                ),

        # depends on StrandPhaseR
#        expand('output/evaluation/variant_phasings/by_approach/{approach1}_vs_{approach2}/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}_{approach1}_vs_{sts_reads}.comp.tsv',
#                approach1=['canonical'],
#                approach2=['strandseq'],
#                var_caller=['longshot', 'freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq']
#                ),

        expand('output/evaluation/known_reference/quastlg_busco/{approach}.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known_ref}.{genemodel}/report.pdf',
                approach=['strandseq'],
                var_caller=['longshot'],
                gq=config['filter_vcf_gq'],
                dp=config['filter_vcf_dp'],
                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
                vc_reads=['HG00733_sra_pbsq1-clr_1000'],
                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                hap=['h1', 'h2', 'h1-un', 'h2-un'],
                known_ref=['GRCh38_1KG_decoy'],
                genemodel=['GRCh38_GENCODEv31_basic']
                ),

        expand('output/evaluation/known_reference/mummer/{approach}.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap}.{known_ref}/{known_ref}_vs_{hap_reads}.{hap}.delta',
                approach=['strandseq'],
                var_caller=['longshot'],
                gq=config['filter_vcf_gq'],
                dp=config['filter_vcf_dp'],
                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
                vc_reads=['HG00733_sra_pbsq1-clr_1000'],
                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
                hap=['h1', 'h2', 'h1-un', 'h2-un'],
                known_ref=['GRCh38_1KG_decoy'],
                ),

        # depends on SaarClust / StrandPhaseR
#        expand('output/evaluation/known_reference/quastlg_busco/{approach}.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known}.{genemodel}/report.pdf',
#                approach=['strandseq'],
#                var_caller=['longshot', 'freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_hgsvc_pbsq2-ccs_clustV1-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                hap=['h1', 'h2', 'h1-un', 'h2-un']
#                ),
