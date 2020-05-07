

rule dump_haploid_read_coverage:
    input:
        'output/alignments/hap_reads_to_reference/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.psort.sam.bam'
    output:
        'output/alignments/hap_reads_to_reference/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.bedGraph'
    log:
        'log/output/alignments/hap_reads_to_reference/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.bg.log'
    benchmark:
        'run/output/alignments/hap_reads_to_reference/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.bg.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt,
        mem_total_mb = 2048,
        mem_per_cpu_mb = 2048
    shell:
        'bedtools genomecov -ibam {input} > {output} 2> {log}'


rule convert_hap_read_coverage:
    input:
       bg_track = 'output/alignments/hap_reads_to_reference/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.bedGraph',
       sizes = 'references/assemblies/{aln_reference}.sizes'
    output:
       'output/cov_tracks/hap_reads/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.bigWig'
    log:
       'log/output/cov_tracks/hap_reads/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.log'
    benchmark:
       'run/output/cov_tracks/hap_reads/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt,
        mem_total_mb = 2048,
        mem_per_cpu_mb = 2048
    shell:
         'bedGraphToBigWig {input.bg_track} {input.sizes} {output} > {log}'
