
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
       expand('output/evaluation/known_reference/quastlg_busco/strandseq.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{assembler}.{pol_reads}.{polisher}.{hap}.{known}.{genemodel}/report.pdf',
               var_caller=['freebayes'],
               assembler=['wtdbg'],
               gq=config['filter_vcf_gq']['freebayes'],
               dp=config['filter_vcf_dp']['freebayes'],
               reference=['HG00733_sra_pbsq1-clr_scV3-wtdbg'],
               vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
               hap_reads=['HG00733_sra_pbsq1-clr_1000'],
               sts_reads=['HG00733_1kg_il25k-npe_sseq'],
               pol_reads=['HG00733_sra_pbsq1-clr_1000'],
               polisher=['arrow-p1'],
               hap=['h1-un', 'h2-un', 'h1', 'h2'],
               known_ref=['GRCh38_GCA_p13'],
               genemodel=['GRCh38_GENCODEv31_basic']
               ),

       expand('output/evaluation/known_reference/quastlg_busco/strandseq.{var_caller}_GQ{gq}_DP{dp}.{reference}.{vc_reads}.{sts_reads}.{hap_reads}.{assembler}.{pol_reads}.{polisher}.{hap}.{known}.{genemodel}/report.pdf',
               var_caller=['freebayes'],
               assembler=['wtdbg'],
               gq=config['filter_vcf_gq']['freebayes'],
               dp=config['filter_vcf_dp']['freebayes'],
               reference=['HG00733_sra_pbsq1-clr_scV3-wtdbg'],
               vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
               hap_reads=['HG00733_sra_pbsq1-clr_1000'],
               sts_reads=['HG00733_1kg_il25k-npe_sseq'],
               pol_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
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
               reference=['HG00733_hgsvc_pbsq2-ccs_scV3-wtdbg'],
               vc_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
               hap_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
               sts_reads=['HG00733_1kg_il25k-npe_sseq'],
               pol_reads=['HG00733_hgsvc_pbsq2-ccs_1000'],
               polisher=['racon-p1'],
               hap=['h1-un', 'h2-un', 'h1', 'h2'],
               known_ref=['GRCh38_GCA_p13'],
               genemodel=['GRCh38_GENCODEv31_basic']
               ),

