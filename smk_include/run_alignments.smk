
localrules: master_run_alignments


rule master_run_alignments:
    input:


rule derive_minimap_parameter_preset:
    input:
        '{filepath}.fastq.gz'
    output:
        '{filepath}.preset.minimap'
    resources:
        runtime_hrs = 0,
        runtime_min = 5
    run:
        preset = ' -x '
        path, filename = os.path.split(input[0])
        tech_spec = filename.split('_')[2]
        if tech_spec.startswith('ont'):
            preset += 'map-ont'
        elif '-clr' in tech_spec:
            # as per recommendation in help
            preset += 'map-pb -H '
        elif '-ccs' in tech_spec:
            # https://github.com/lh3/minimap2/issues/325
            preset += 'asm20'
        else:
            raise ValueError('No minimap preset for file: {}'.format(filename))

        # This is special for alignment before Racon polishing;
        # note that even for CCS reads, the PacBio preset is used
        # Recommended by Aaron Wenger
        if all([x in path for x in ['diploid_assembly', 'draft', 'haploid_fastq']]):
            preset = ' -x map-pb --eqx -m 5000 --secondary=no '

        with open(output[0], 'w') as dump:
            _ = dump.write(preset)


rule derive_pbmm2_parameter_preset:
    input:
        '{filepath}.pbn.bam'
    output:
        '{filepath}.preset.pbmm2'
    resources:
        runtime_hrs = 0,
        runtime_min = 5
    run:
        filename = os.path.basename(input[0])
        tech_spec = filename.split('_')[2]
        preset = ''
        if '-clr' in tech_spec:
            preset = 'SUBREAD'
        elif '-ccs' in tech_spec:
            preset = 'CCS'
        else:
            raise ValueError('No pbmm2 preset for file: {}'.format(filename))

        with open(output[0], 'w') as dump:
            _ = dump.write(preset)


rule minimap_reads_to_reference_alignment:
    input:
        reads = 'input/fastq/complete/{sample}.fastq.gz',
        preset = 'input/fastq/complete/{sample}.preset.minimap',
        reference = 'output/reference_assembly/{folder_path}/{reference}.fasta'
    output:
        'output/alignments/reads_to_reference/{folder_path}/{sample}_map-to_{reference}.sam'
    log:
        'log/output/alignments/reads_to_reference/{folder_path}/{sample}_map-to_{reference}.log'
    benchmark:
        'run/output/alignments/reads_to_reference/{{folder_path}}/{{sample}}_map-to_{{reference}}.t{}.rsrc'.format(config['num_cpu_max'])
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = int(49152 / config['num_cpu_max']),
        mem_total_mb = 49152,
        runtime_hrs = 3
    params:
        individual = lambda wildcards: wildcards.sample.split('_')[0],
        preset = load_preset_file
    shell:
        'minimap2 -t {threads} {params.preset} -a -o {output} ' \
            '-R "@RG\\tID:1\\tSM:{params.individual}" ' \
            '{input.reference} {input.reads} &> {log}'


rule pbmm2_reads_to_reference_alignment:
    input:
        reads = 'input/bam/complete/{sample}.pbn.bam',
        preset = 'input/bam/complete/{sample}.preset.pbmm2',
        reference = 'output/reference_assembly/{folder_path}/{reference}.fasta'
    output:
        bam = 'output/alignments/reads_to_reference/{folder_path}/{sample}_map-to_{reference}.psort.pbn.bam',
    log:
        'log/output/alignments/reads_to_reference/{folder_path}/{sample}_map-to_{reference}.pbn.log'
    benchmark:
        'run/output/alignments/reads_to_reference/{{folder_path}}/{{sample}}_map-to_{{reference}}.pbn.t{}.rsrc'.format(config['num_cpu_max'])
    conda:
        config['conda_env_pbtools']
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = int(98304 / config['num_cpu_max']),
        mem_total_mb = 98304,
        runtime_hrs = 23
    params:
        align_threads = config['num_cpu_max'] - config['num_cpu_low'],
        sort_threads = config['num_cpu_low'],
        sort_memory_mb = int((98304 / config['num_cpu_max']) * config['num_cpu_low']),
        preset = load_preset_file,
        individual = lambda wildcards: wildcards.sample.split('_')[0]
    shell:
        'pbmm2 align --log-level INFO --sort --sort-memory {params.sort_memory_mb}M --no-bai ' \
            ' --alignment-threads {params.align_threads} --sort-threads {params.sort_threads} ' \
            ' --preset {params.preset} --sample {params.individual} ' \
            ' {input.reference} {input.reads} {output.bam} &> {log}'


def select_bwa_index(wildcards):

    sts_individual = wildcards.sts_reads.split('_')[0]
    ref_individual = wildcards.reference.split('_')[0]

    assert sts_individual == ref_individual, 'Mixed individual match (bwa index): {}'.format(wildcards)

    if '_nhr-' in wildcards.reference:
        # non-haplotype resolved assembly / collapsed assembly
        idx = 'output/reference_assembly/non-hap-res/{}/bwa_index/{}.bwt'.format(wildcards.reference, wildcards.reference)
    elif '_scV{}-'.format(config['git_commit_version']) in wildcards.reference:
        idx = 'output/reference_assembly/clustered/{}/{}/bwa_index/{}.bwt'.format(wildcards.sts_reads, wildcards.reference, wildcards.reference)
    else:
        raise ValueError('Unexpected reference type: {} / {}'.format(wildcards.reference, wildcards))
    return idx


rule bwa_strandseq_to_reference_alignment:
    input:
        mate1 = 'input/fastq/strand-seq/{sts_reads}/{individual}_{sample_id}_1.fastq.gz',
        mate2 = 'input/fastq/strand-seq/{sts_reads}/{individual}_{sample_id}_2.fastq.gz',
        ref_index = select_bwa_index,
        sts_reads = 'input/fastq/strand-seq/{sts_reads}.fofn'
    output:
        bam = 'output/alignments/strandseq_to_reference/{reference}/{sts_reads}/temp/aln/{individual}_{sample_id}.filt.sam.bam'
    log:
        bwa = 'log/output/alignments/strandseq_to_reference/{reference}/{sts_reads}/{individual}_{sample_id}.bwa.log',
        samtools = 'log/output/alignments/strandseq_to_reference/{reference}/{sts_reads}/{individual}_{sample_id}.samtools.log',
    benchmark:
        'run/output/alignments/strandseq_to_reference/{{reference}}/{{sts_reads}}/{{individual}}_{{sample_id}}.t{}.rsrc'.format(config['num_cpu_low'])
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = int(8192 / config['num_cpu_low']),
        mem_total_mb = 8192
    params:
        idx_prefix = lambda wildcards, input: input.ref_index.split('.')[0]
    shell:
        'bwa mem -t {threads}' \
            ' -R "@RG\\tID:{wildcards.individual}_{wildcards.sample_id}\\tPL:Illumina\\tSM:{wildcards.individual}"' \
            ' -v 2 {params.idx_prefix} {input.mate1} {input.mate2} 2> {log.bwa} | ' \
            ' samtools view -b -F 2304 /dev/stdin > {output.bam} 2> {log.samtools}'


rule minimap_racon_polish_alignment_pass1:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    pol_reads = FASTQ file used for Racon contig polishing
    """
    input:
        reads = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{pol_reads}.{hap}.{sequence}.fastq.gz',
        preset = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{pol_reads}.{hap}.{sequence}.preset.minimap',
        contigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.fasta',
    output:
        sam = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.racon-p1.psort.sam',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.racon-p1.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT_PROTECTED + '/polishing/alignments/{{pol_reads}}_map-to_{{hap_reads}}-{{assembler}}.{{hap}}.{{sequence}}.racon-p1.t{}.rsrc'.format(config['num_cpu_high'])
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(attempt * 16384 / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: attempt * 16384,
        mem_sort_mb = 4096,
        runtime_hrs = 3
    params:
        individual = lambda wildcards: wildcards.hap_reads.split('_')[0],
        preset = load_preset_file,
        discard_flag = config['minimap_racon_aln_discard'],
        min_qual = config['minimap_racon_aln_min_qual']
    shell:
        'minimap2 -t {threads} -a {params.preset} -R "@RG\\tID:1\\tSM:{params.individual}" ' \
            ' {input.contigs} {input.reads} 2> {log} | samtools sort -m {resources.mem_sort_mb}M | ' \
            ' samtools view -q {params.min_qual} -F {params.discard_flag} /dev/stdin > {output.sam}'


rule pbmm2_arrow_polish_alignment_pass1:
    input:
        reads = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pol_reads}.{hap}.{sequence}.pbn.bam',
        preset = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pol_reads}.{hap}.{sequence}.preset.pbmm2',
        contigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.fasta',
    output:
        bam = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.psort.pbn.bam',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.psort.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT_PROTECTED + '/polishing/alignments/{{pol_reads}}_map-to_{{hap_reads}}-{{assembler}}.{{hap}}.{{sequence}}.arrow-p1.psort.t{}.rsrc'.format(config['num_cpu_high'])
    conda:
        config['conda_env_pbtools']
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(attempt * 16384 / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: attempt * 16384,
        runtime_hrs = 4
    params:
        align_threads = config['num_cpu_high'] - 2,
        sort_threads = 2,
        preset = load_preset_file,
        individual = lambda wildcards: wildcards.hap_reads.split('_')[0]
    shell:
        'pbmm2 align --log-level INFO --sort --sort-memory {resources.mem_per_cpu_mb}M --no-bai ' \
            ' --alignment-threads {params.align_threads} --sort-threads {params.sort_threads} ' \
            ' --preset {params.preset} --sample {params.individual} ' \
            ' {input.contigs} {input.reads} {output.bam} &> {log}'
