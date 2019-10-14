
include: 'eval_known_reference.smk'
include: 'eval_variant_calls.smk'
include: 'arrow_polishing.smk'

localrules: master_results_child

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
       expand('output/evaluation/known_reference/quastlg_busco/strandseq.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{pol_reads}.{polisher}.{hap}.{known_ref}.{genemodel}/report.pdf',
               var_caller=['freebayes'],
               gq=config['filter_vcf_gq']['freebayes'],
               dp=config['filter_vcf_dp']['freebayes'],
               reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
               vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
               hap_reads=['HG00733_sra_pbsq1-clr_1000'],
               sts_reads=['HG00733_1kg_il25k-npe_sseq'],
               pol_reads=['HG00733_sra_pbsq1-clr_1000'],
               polisher=['arrow-p1'],
               hap=['h1-un', 'h2-un'],
               known_ref=['GRCh38_GCA_p13'],
               genemodel=['GRCh38_GENCODEv31_basic']
               ),

#        # Running on UMBRIEL
#        expand('output/evaluation/known_reference/quastlg_busco/strandseq.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{pol_reads}.{polisher}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                pol_reads=['HG00733_sra_pbsq1-clr_1000'],
#                polisher=['arrow-p1'],
#                hap=['h2-un'],
#                known_ref=['GRCh38_GCA_p13'],
#                genemodel=['GRCh38_GENCODEv31_basic']
#                ),

#        # Running on CORDELIA
#        expand('output/evaluation/known_reference/quastlg_busco/strandseq.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{pol_reads}.{polisher}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                pol_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                polisher=['racon-p1'],
#                hap=['h1-un', 'h2-un'],
#                known_ref=['GRCh38_GCA_p13'],
#                genemodel=['GRCh38_GENCODEv31_basic']
#                ),

#        # Running on MIRANDA
#        expand('output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{reference}.{hap}.{polisher}.psort.pbn.bam',
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                pol_reads=['HG00733_sra_pbsq1-clr_1000'],
#                polisher=['arrow-p1'],
#                hap=['h2-un'],
#                ),

#        # Running on BIANCA
#        expand('output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{reference}.{hap}.{polisher}.psort.pbn.bam',
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                pol_reads=['HG00733_sra_pbsq1-clr_1000'],
#                polisher=['arrow-p1'],
#                hap=['h1-un'],
#                ),

#        # TODO on CORDELIA
#        expand('output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/{hap_reads}.{hap}.{polisher}.fasta',
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV2-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                pol_reads=['HG00733_sra_pbsq1-clr_1000'],
#                polisher=['arrow-p1'],
#                hap=['h2-un']
#                ),

#        # Running on CORDELIA
#        expand('output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/{hap_reads}.{hap}.arrow-p1.{sequence}.fasta',
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV2-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                pol_reads=['HG00733_sra_pbsq1-clr_1000'],
#                polisher=['arrow-p1'],
#                hap=['h2-un'],
#                sequence=list(['ctg' + str(i) for i in range(1, 3354)])
#                ),

#        # ERROR on UMBRIEL
#        expand('output/evaluation/known_reference/quastlg_busco/strandseq.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{pol_reads}.{polisher}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV2-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                pol_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                polisher=['racon-p1'],
#                hap=['h1', 'h2', 'h1-un', 'h2-un'],
#                known_ref=['GRCh38_GCA_p13'],
#                genemodel=['GRCh38_GENCODEv31_basic']
#                ),

#        # TODO running on JULIET
#        expand('output/evaluation/known_reference/quastlg_busco/strandseq.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{pol_reads}.{polisher}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                pol_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                polisher=['racon-p1'],
#                hap=['h1', 'h2', 'h1-un', 'h2-un'],
#                known_ref=['GRCh38_GCA_p13'],
#                genemodel=['GRCh38_GENCODEv31_basic']
#                ),

#        # DONE
#        expand('output/evaluation/known_reference/quastlg_busco/strandseq.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{pol_reads}.{polisher}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_hgsvc_pbsq2-ccs_clustV2-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                pol_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                polisher=['racon-p1'],
#                hap=['h1', 'h2', 'h1-un', 'h2-un'],
#                known_ref=['GRCh38_GCA_p13'],
#                genemodel=['GRCh38_GENCODEv31_basic']
#                ),

#        # DONE
#        expand('output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{reference}.{hap}.arrow-p1.psort.pbn.bam',
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                pol_reads=['HG00733_sra_pbsq1-clr_1000'],
#                hap=['h1', 'h2', 'h1-un', 'h2-un'],
#                ),

#        # DONE
#        expand('output/alignments/reads_to_reference/{sample}_map-to_{reference}.psort.pbn.bam',
#                sample=['HG00733_sra_pbsq1-clr_1000'],
#                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb', 'HG00733_sra_pbsq1-clr_clustV2-100kb']),

#          # DONE
#         expand('output/evaluation/known_reference/quastlg_busco/{approach}.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                 approach=['strandseq'],
#                 var_caller=['freebayes'],
#                 gq=config['filter_vcf_gq'],
#                 dp=config['filter_vcf_dp'],
#                 reference=['HG00733_hgsvc_pbsq2-ccs_clustV2-100kb'],
#                 vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                 hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                 sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                 hap=['h1', 'h2', 'h1-un', 'h2-un'],
#                 known_ref=['GRCh38_GCA_p13'],
#                 genemodel=['GRCh38_GENCODEv31_basic']
#                 ),

#        # DONE
#        expand('output/evaluation/known_reference/quastlg_busco/{approach}.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                approach=['strandseq'],
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV2-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                hap=['h1', 'h2', 'h1-un', 'h2-un'],
#                known_ref=['GRCh38_GCA_p13'],
#                genemodel=['GRCh38_GENCODEv31_basic']
#                ),

#        # DONE
#        expand(rules.merge_sequence_vcf_files.output,
#                var_caller=['freebayes'],
#                reference=['HG00733_sra_pbsq1-clr_clustV2-100kb'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'])

#        # DONE
#        expand('output/saarclust/results/{individual}_hgsvc_pbsq2-ccs_sqa/{individual}_1kg_il25k-npe_sseg',
#                individual=['HG00733'])

#         # DONE
#         'output/variant_calls/longshot/HG00733_sra_pbsq1-clr_clustV2-100kb/final_GQ100_DP50/HG00733_sra_pbsq1-clr_1000.final.vcf'

#        # DONE
#        expand('output/evaluation/known_reference/quastlg_busco/{approach}.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{hap_reads}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                approach=['canonical'],
#                var_caller=['longshot'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['hg38_GCA_p13'],
#                vc_reads=['HG00733_sra_pbsq1-clr_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                hap=['h1', 'h2', 'h1-un', 'h2-un'],
#                known_ref=['GRCh38_GCA_p13'],
#                genemodel=['GRCh38_GENCODEv31_basic']
#                ),

#        # DONE
#        expand('output/saarclust/results/{individual}_hgsvc_pbsq2-ccs_sqa/{individual}_1kg_il25k-npe_sseg',
#                individual=['HG00733', 'HG00732', 'HG00731'])

#        # DONE
#        expand('output/evaluation/known_reference/quastlg_busco/{approach}.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                approach=['strandseq'],
#                var_caller=['freebayes'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
#                vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                hap=['h1', 'h2', 'h1-un', 'h2-un'],
#                known_ref=['GRCh38_GCA_p13'],
#                genemodel=['GRCh38_GENCODEv31_basic']
#                ),

#        # DONE
#        expand('output/evaluation/known_reference/quastlg_busco/{approach}.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{hap}.{known_ref}.{genemodel}/report.pdf',
#                approach=['strandseq'],
#                var_caller=['longshot'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
#                vc_reads=['HG00733_sra_pbsq1-clr_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                hap=['h1', 'h2', 'h1-un', 'h2-un'],
#                known_ref=['GRCh38_GCA_p13'],
#                genemodel=['GRCh38_GENCODEv31_basic']
#                ),
#
#        # DONE
#        expand('output/evaluation/known_reference/mummer/{approach}.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap}.{known_ref}/{known_ref}_vs_{hap_reads}.{hap}.delta',
#                approach=['strandseq'],
#                var_caller=['longshot'],
#                gq=config['filter_vcf_gq'],
#                dp=config['filter_vcf_dp'],
#                reference=['HG00733_sra_pbsq1-clr_clustV1-100kb'],
#                vc_reads=['HG00733_sra_pbsq1-clr_1000'],
#                hap_reads=['HG00733_sra_pbsq1-clr_1000'],
#                sts_reads=['HG00733_1kg_il25k-npe_sseq'],
#                hap=['h1', 'h2', 'h1-un', 'h2-un'],
#                known_ref=['GRCh38_GCA_p13'],
#                ),
