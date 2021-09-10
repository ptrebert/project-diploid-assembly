
def set_alignment_memory(wildcards, attempt):
    """
    aligning the ONT-UL reads is driving me crazy...
    now use the hammer
    """
    base_mem = 94208
    if 'ONTEC' in wildcards.read_type:
        base_mem += 24768
    elif 'ONTUL' in wildcards.read_type:
        base_mem = 262144
    elif 'HIFIEC' in wildcards.read_type:
        pass
    elif 'HIFIAF' in wildcards.read_type:
        pass
    else:
        raise ValueError(str(wildcards))
    return base_mem * attempt


def set_alignment_runtime(wildcards, attempt):
    """
    aligning the ONT-UL reads is driving me crazy...
    now use the hammer
    """
    base_hrs = 36
    if 'ONTEC' in wildcards.read_type:
        base_hrs = 72
    elif 'ONTUL' in wildcards.read_type:
        base_hrs = 167
    elif 'HIFIEC' in wildcards.read_type:
        pass
    elif 'HIFIAF' in wildcards.read_type:
        pass
    else:
        raise ValueError(str(wildcards))
    runtime_limit = base_hrs * attempt
    if runtime_limit > 167:
        raise ValueError(f'Runtime limit exceeds longq limit: {runtime_limit} / {str(wildcards)}')
    return runtime_limit


def set_alignment_preset(wildcards):
    preset = None
    if 'ONTEC' in wildcards.read_type:
        preset = 'map-pb'
    elif 'ONTUL' in wildcards.read_type:
        preset = 'map-ont'
    elif 'HIFIEC' in wildcards.read_type:
        preset = 'map-pb'
    elif 'HIFIAF' in wildcards.read_type:
        preset = 'map-pb'
    else:
        raise ValueError(str(wildcards))
    assert preset is not None
    return preset


rule qc_mmap_align_readsets:
    """
    https://github.com/lh3/minimap2/issues/771
    Above github issue contains some hints how to speed up alignment
    for ONT reads. Following this, set...
    -k17
    --cap-kalloc=1g
    -K4g
    """
    input:
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.read_type],
        reference = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta'),
    output:
        bam = 'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.psort.bam',
        bai = 'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.psort.bam.bai'
    log:
        'log/output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.psort.mmap.log'
    benchmark:
        'rsrc/output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.psort.mmap.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: set_alignment_memory(wildcards, attempt),
        runtime_hrs = lambda wildcards, attempt: set_alignment_runtime(wildcards, attempt),
        mem_sort_mb = 4096,
        align_threads = config['num_cpu_high'] - config['num_cpu_low'],
        sort_threads = config['num_cpu_low'],
    params:
        individual = lambda wildcards: wildcards.sample,
        readgroup_id = lambda wildcards: f'{wildcards.read_type}_{wildcards.readset.replace('.', '')}',
        preset = lambda wildcards: set_alignment_preset(wildcards),
        temp_prefix = lambda wildcards: f'temp/mmap/{wildcards.reference}/{wildcards.sample}/{wildcards.read_type}/{wildcards.readset}/tmp_stsort_',
        temp_dir = lambda wildcards: f'temp/mmap/{wildcards.reference}/{wildcards.sample}/{wildcards.read_type}/{wildcards.readset}/',
        validate = lambda wildcards, input:validate_readset(wildcards.readset, input.reads)
    shell:
        'rm -rfd {params.temp_dir} && mkdir -p {params.temp_dir} && '
        'minimap2 -t {resources.align_threads} -Y -L --eqx --MD -a -x {params.preset} '
        '-R "@RG\\tID:{params.readgroup_id}\\tSM:{params.individual}" --secondary=no '
        '-k17 --cap-kalloc=1g -K4g '
        '{input.reference} {input.reads} | '
        'samtools view -u -F 260 | '
        'samtools sort -m {resources.mem_sort_mb}M -@ {resources.sort_threads} -T {params.temp_prefix} -O BAM > {output.bam} '
        ' && '
        'samtools index {output.bam} ; rm -rfd {params.temp_dir}'


rule dump_alignment_coverage:
    input:
        bam = 'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.psort.bam',
        bai = 'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.psort.bam.bai',
    output:
        bed = 'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.aln.bed',
    benchmark:
        'rsrc/output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.aln.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 1
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt ** 3,
    shell:
        'bedtools bamtobed -i {input.bam} > {output.bed}'
