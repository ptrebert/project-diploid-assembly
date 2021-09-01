
rule qc_merge_ontul:
    input:
        fastq = lambda wildcards: SAMPLE_INFOS[wildcards.sample]['ONTUL']
    output:
        'input/ONTUL/{sample}_ONTUL_guppy-5.0.11-sup-prom.fasta.gz'
    benchmark:
        'rsrc/input/ONTUL/{sample}_ONTUL_guppy-5.0.11-sup-prom.merge.rsrc'
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    shell:
        'pigz -d -c {input} | seqtk seq -A -C | pigz -p {threads} --best > {output}'


rule qc_local_seq_stats:
    input:
        'input/{read_type}/{sample}_{read_type}_{readset}.fasta.gz'
    output:
        'input/{read_type}/{sample}_{read_type}_{readset}.stats.tsv.gz'
    conda: '../../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        read_type = '(ONTUL|ONTEC|ONTHY)'
    threads: 2
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    shell:
        'pigz -d -c {input} | seqtk comp | pigz --best > {output}'


rule qc_remote_seq_stats:
    input:
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.read_type]
    output:
        'input/{read_type}/{sample}_{read_type}_{readset}.stats.tsv.gz'
    conda: '../../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        read_type = '(HIFIEC|HIFIAF|SHORT)'
    threads: 2
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    shell:
        'pigz -d -c {input} | seqtk comp | pigz --best > {output}'
