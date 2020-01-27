
# just to get rid off 'unresolved reference' indicators
import os


rule racon_contig_polishing_pass1:
    input:
        pol_reads = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{pol_reads}.{hap}.{sequence}.fastq.gz',
        contigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.fasta',
        seq_info = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.fasta.fai',
        alignments = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.racon-p1.psort.sam',
    output:
        'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.racon-p1.fasta'
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.racon-p1.log'
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.racon-p1.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((12288 * attempt) / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 2 * attempt
    shell:
        'racon --threads {threads} --include-unpolished {input.pol_reads} {input.alignments} {input.contigs} > {output} 2> {log}'


rule racon_contig_polishing_pass2:
    input:
        pol_reads = os.path.join('output', PATH_STRANDSEQ_DGA_SPLIT, 'draft/haploid_fastq/{pol_reads}.{hap}.{sequence}.fastq.gz'),
        contigs = os.path.join('output', PATH_STRANDSEQ_DGA_SPLIT, 'polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.racon-p1.fasta'),
        seq_info = os.path.join('output', PATH_STRANDSEQ_DGA_SPLIT, 'polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.racon-p1.fasta.fai'),
        alignments = os.path.join('output', PATH_STRANDSEQ_DGA_SPLIT, 'polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.racon-p2.psort.sam'),
    output:
        os.path.join('output', PATH_STRANDSEQ_DGA_SPLIT, 'polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.racon-p2.fasta')
    log:
        os.path.join('log/output', PATH_STRANDSEQ_DGA_SPLIT, 'polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.racon-p2.log')
    benchmark:
        os.path.join('run/output', PATH_STRANDSEQ_DGA_SPLIT, 'polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.racon-p2.' + 't{}.rsrc'.format(config['num_cpu_medium']))
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((12288 * attempt) / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 2 * attempt
    shell:
        'racon --threads {threads} --include-unpolished {input.pol_reads} {input.alignments} {input.contigs} > {output} 2> {log}'
