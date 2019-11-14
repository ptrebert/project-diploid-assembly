
include: 'canonical_dga.smk'
include: 'strandseq_dga_joint.smk'
include: 'strandseq_dga_split.smk'
include: 'aux_utilities.smk'
include: 'run_alignments.smk'


rule racon_contig_polishing_pass1:
    input:
        pol_reads = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{pol_reads}.{hap}.{sequence}.fastq.gz',
        contigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}.{sequence}.fasta',
        seq_info = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}.{sequence}.fasta.fai',
        alignments = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.psort.sam',
    output:
        'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{sequence}.racon-p1.fasta'
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{sequence}.racon-p1.log'
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{sequence}.racon-p1.rsrc'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = int(24576 / config['num_cpu_medium']),
        mem_total_mb = 24576,
        runtime_hrs = 8
    shell:
        'racon --threads {threads} --include-unpolished {input.pol_reads} {input.alignments} {input.contigs} > {output} 2> {log}'
