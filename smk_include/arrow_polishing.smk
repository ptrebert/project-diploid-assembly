
# rule arrow_contig_polishing_pass1:
#     input:
#         contigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.fasta',
#         seq_info = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.fasta.fai',
#         alignments = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.psort.pbn.bam',
#         aln_index = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.psort.pbn.bam.pbi',
#     output:
#         'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.fasta'
#     log:
#         'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.log'
#     benchmark:
#         'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.rsrc'
#     conda:
#          '../environment/conda/conda_pbtools.yml'
#     threads: config['num_cpu_medium']
#     resources:
#         mem_per_cpu_mb = lambda wildcards, attempt: int(12288 + attempt * 24576 / config['num_cpu_medium']),
#         mem_total_mb = lambda wildcards, attempt: 12288 + attempt * 24576,
#         runtime_hrs = 47
#     shell:
#         'variantCaller --algorithm=arrow --log-file {log} --log-level INFO -j {threads} '
#             ' --reference {input.contigs} -o {output} {input.alignments}'


rule arrow_contig_polishing_pass1:
    input:
         contigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.fasta',
         seq_info = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.fasta.fai',
         alignments = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.psort.pbn.bam',
         aln_index = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.psort.pbn.bam.bai',
    output:
          fasta = os.path.join('output', PATH_STRANDSEQ_DGA_SPLIT, 'polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.fasta'),
          gff = os.path.join('output', PATH_STRANDSEQ_DGA_SPLIT, 'polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.gff'),
    log:
       'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.log'
    benchmark:
        os.path.join('run/output', PATH_STRANDSEQ_DGA_SPLIT, 'polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1' + '.t{}.rsrc'.format(config['num_cpu_medium']))
    conda:
         '../environment/conda/conda_pbtools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(12288 + attempt * 24576 / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 12288 + attempt * 24576,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    shell:
         'gcpp --num-threads {threads} --algorithm=arrow  --log-level INFO --log-file {log} '
         '--reference {input.contigs} --output {output.fasta},{output.gff} {input.alignments}'
