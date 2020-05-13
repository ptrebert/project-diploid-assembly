

rule dump_haploid_read_coverage:
    """
    "Since recently", UCSC tools require old "ASCII" sort order for the big* indices
    to be correct. This is incompatible with default locale (UTF-8) on many Linux systems
    """
    input:
        'output/alignments/hap_reads_to_reference/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.psort.sam.bam'
    output:
        'output/alignments/hap_reads_to_reference/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.sorted.bedGraph'
    log:
        bedtools = 'log/output/alignments/hap_reads_to_reference/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.bg.log',
        sort = 'log/output/alignments/hap_reads_to_reference/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.sort.log',
    benchmark:
        os.path.join('run/output/alignments/hap_reads_to_reference',
                     '{folder_path}',
                     '{file_name}_map-to_{aln_reference}.{hap}.bg' + '.t{}.rsrc'.format(config['num_cpu_medium']))
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * 6,
        mem_total_mb = lambda wildcards, attempt: attempt * 32768 + 32768,
        mem_per_cpu_mb = lambda wildcards, attempt: int((attempt * 32768 + 32768) / config['num_cpu_medium'])
    shell:
        'bedtools genomecov -bg -ibam {input} 2> {log.bedtools}'
        ' | '
        'LC_COLLATE=C sort --buffer-size={resources.mem_total_mb}M --parallel={threads} '
        '-k1,1 -k2,2n > {output} 2> {log.sort}'


rule convert_hap_read_coverage:
    input:
       bg_track = 'output/alignments/hap_reads_to_reference/{folder_path}/{file_name}_map-to_{aln_reference}.{hap}.sorted.bedGraph',
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
        runtime_hrs = lambda wildcards, attempt: attempt * 6,
        mem_total_mb = lambda wildcards, attempt: attempt * 16384,
        mem_per_cpu_mb = lambda wildcards, attempt: attempt * 16384
    shell:
         'bedGraphToBigWig {input.bg_track} {input.sizes} {output} 2> {log}'
