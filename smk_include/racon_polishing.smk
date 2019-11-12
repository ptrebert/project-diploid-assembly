
include: 'canonical_dga.smk'
include: 'strandseq_dga_joint.smk'
include: 'strandseq_dga_split.smk'
include: 'run_alignments.smk'


rule racon_contig_polishing_pass1:
    input:
        pol_reads = 'output/diploid_assembly/{folder_path}/draft/haploid_fastq/{pol_reads}.{hap}{split}.fastq.gz',
        alignments = 'output/diploid_assembly/{folder_path}/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}{split}.racon-p1.psort.sam',
        contigs = 'output/diploid_assembly/{folder_path}/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}{split}.fasta'
    output:
        'output/diploid_assembly/{folder_path}/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}{split}.racon-p1.fasta'
    log:
        'log/output/diploid_assembly/{folder_path}/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}{split}.racon-p1.log'
    benchmark:
        'run/output/diploid_assembly/{folder_path}/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}{split}.racon-p1.rsrc'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = 14336,
        mem_total_mb = 172032
    shell:
        'racon --threads {threads} --include-unpolished {input.pol_reads} {input.alignments} {input.contigs} > {output} 2> {log}'
