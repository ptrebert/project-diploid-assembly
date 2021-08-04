include: 'assmy_collect_input.smk'


rule dump_tigs_to_fasta:
    input:
        tigs = lambda wildcards: str(SAMPLE_INFOS[wildcards.sample][wildcards.tigs]),
    output:
        fasta = 'output/tigs/{sample_info}_{sample}.{tigs}.fasta'
    benchmark:
        'output/tigs/{sample_info}_{sample}.{tigs}.dump.rsrc'
    wildcard_constraints:
        tigs = '(TIGPRI|TIGALT|TIGRAW)',
        sample = CONSTRAINT_SAMPLES
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        header_prefix = f'>{wildcards.sample_info}_{wildcards.sample}.'
        with open(output.fasta, 'w') as fasta:
            with open(input.tigs, 'r') as gfa:
                for ln, line in enumerate(gfa, start=1):
                    if not line.startswith('S'):
                        continue
                    try:
                        parts = line.split()
                        node_id = parts[1].strip()
                        sequence = parts[2].strip()
                    except (ValueError, IndexError) as err:
                        if not line.strip():
                            continue
                        raise ValueError(f'LN {ln} : FILE {gfa_file}: {str(err)}')
                    _ = fasta.write(f'{header_prefix}{node_id}\n{sequence}\n')
    # END OF RUN BLOCK


rule align_tigs_to_reference:
    input:
        fa_tigs = 'output/tigs/{sample_info}_{sample}.{tigs}.fasta',
        fa_ref = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta')
    output:
        bam = 'output/tig_aln/split/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}.psort.bam',
        bai = 'output/tig_aln/split/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}.psort.bam.bai',
        temp_dir = temp(
            directory('temp/minimap/tig_aln/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}')
        )
    log:
        mm = 'log/output/tig_aln/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}.mm.log',
        view = 'log/output/tig_aln/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}.view.log',
        sort = 'log/output/tig_aln/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}.sort.log',
    benchmark:
        'rsrc/output/tig_aln/{sample_info}_{sample}_{tigs}_MAP-TO_{reference}.rsrc',
    wildcard_constraints:
        tigs = '(TIGPRI|TIGALT|TIGRAW)'
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: (config['num_cpu_medium'] + config['num_cpu_low'])
    resources:
        mem_total_mb = lambda wildcards, attempt: 36864 + 36864 * attempt,
        runtime_hrs = lambda wildcards, attempt: 16 * attempt,
        mem_sort_mb = 4096
    params:
        align_threads = config['num_cpu_medium'],
        sort_threads = config['num_cpu_low'],
        individual = lambda wildcards: wildcards.sample,
        readgroup_id = lambda wildcards: wildcards.tigs,
        discard_flag = 260  # unmapped or not primary
    shell:
        'mkdir -p {output.temp_dir} && '
        'minimap2 -t {params.align_threads} '
            '--secondary=no --eqx -L -Y -ax asm20 -m 10000 '
            '-R "@RG\\tID:{params.readgroup_id}\\tSM:{params.individual}\\tCN:HHU" '
            '{input.fa_ref} {input.fa_tigs} 2> {log.mm} | '
        'samtools view -u -F {params.discard_flag} /dev/stdin 2> {log.view} | '
        'samtools sort -@ {params.sort_threads} -m {resources.mem_sort_mb}M -T {output.temp_dir}/tmp_part -O BAM -l 6 /dev/stdin > {output.bam} 2> {log.sort} '
        ' && '
        'samtools index {output.bam}'


rule merge_tigs_to_reference_bams:
    input:
        bams = expand(
            'output/tig_aln/split/{{sample_info}}_{{sample}}_{tigs}_MAP-TO_{{reference}}.psort.bam',
            tigs=['TIGPRI', 'TIGALT']
        ),
        bais = expand(
            'output/tig_aln/split/{{sample_info}}_{{sample}}_{tigs}_MAP-TO_{{reference}}.psort.bam.bai',
            tigs=['TIGPRI', 'TIGALT']
        )
    output:
        bam = 'output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.psort.bam',
        bai = 'output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.psort.bam.bai',
    wildcard_constraints:
        sample = CONSTRAINT_SAMPLES
    log:
        'log/output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.merge.log'
    benchmark:
        'rsrc/output/tig_aln/{sample_info}_{sample}_MAP-TO_{reference}.merge.rsrc'
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'samtools merge -l 6 -O BAM --no-PG -@ {threads} {output.bam} {input.bams} 2> {log} '
            ' && '
            'samtools index {output.bam}'