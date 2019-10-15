
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
       expand('output/evaluation/known_reference/quastlg_busco/strandseq.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{assembler}.{pol_reads}.{polisher}.{hap}.{known}.{genemodel}/report.pdf',
               var_caller=['freebayes'],
               assembler=['wtdbg'],
               gq=config['filter_vcf_gq']['freebayes'],
               dp=config['filter_vcf_dp']['freebayes'],
               reference=['HG00731_hgsvc_pbsq2-ccs_scV3-wtdbg'],
               vc_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
               hap_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
               sts_reads=['HG00731_1kg_il25k-npe_sseq'],
               pol_reads=['HG00731_hgsvc_pbsq2-ccs_1000'],
               polisher=['racon-p1'],
               hap=['h1-un', 'h2-un', 'h1', 'h2'],
               known_ref=['GRCh38_GCA_p13'],
               genemodel=['GRCh38_GENCODEv31_basic']
               ),

       expand('output/evaluation/known_reference/quastlg_busco/strandseq.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{assembler}.{pol_reads}.{polisher}.{hap}.{known}.{genemodel}/report.pdf',
               var_caller=['freebayes'],
               assembler=['wtdbg'],
               gq=config['filter_vcf_gq']['freebayes'],
               dp=config['filter_vcf_dp']['freebayes'],
               reference=['HG00732_hgsvc_pbsq2-ccs_scV3-wtdbg'],
               vc_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
               hap_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
               sts_reads=['HG00732_1kg_il25k-npe_sseq'],
               pol_reads=['HG00732_hgsvc_pbsq2-ccs_1000'],
               polisher=['racon-p1'],
               hap=['h1-un', 'h2-un', 'h1', 'h2'],
               known_ref=['GRCh38_GCA_p13'],
               genemodel=['GRCh38_GENCODEv31_basic']
               ),
