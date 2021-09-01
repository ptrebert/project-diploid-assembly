
def select_winnowmap_reads(wildcards):

    if 'MAP-TO_HIFIEC' in wildcards.readset:
        reads = [f'input/ONTEC/{wildcards.sample}_{wildcards.readset}.fasta.gz']
    elif 'ONTUL' in wildcards.readset:
        reads = [f'input/ONTUL/{wildcards.sample}_{wildcards.readset}.fasta.gz']
    elif 'HIFIEC' in wildcards.readset:
        reads = SAMPLE_INFOS[wildcards.sample]['HIFIEC']
    elif 'HIFIAF' in wildcards.readset:
        reads = SAMPLE_INFOS[wildcards.sample]['HIFIAF']
    else:
        raise ValueError(str(wildcards))
    return reads


def set_winnowmap_memory(wildcards, attempt):
    """
    aligning the ONT-UL reads is driving me crazy...
    now use the hammer
    """
    base_mem = 65536
    if 'MAP-TO_HIFIEC' in wildcards.readset:
        base_mem += 24768
    elif 'ONTUL' in wildcards.readset:
        base_mem = 262144
    elif 'HIFIEC' in wildcards.readset:
        pass
    elif 'HIFIAF' in wildcards.readset:
        pass
    else:
        raise ValueError(str(wildcards))
    return base_mem * attempt


def set_winnowmap_runtime(wildcards, attempt):
    """
    aligning the ONT-UL reads is driving me crazy...
    now use the hammer
    """
    base_hrs = 24
    if 'MAP-TO_HIFIEC' in wildcards.readset:
        base_hrs = 60
    elif 'ONTUL' in wildcards.readset:
        base_hrs = 84
    elif 'HIFIEC' in wildcards.readset:
        pass
    elif 'HIFIAF' in wildcards.readset:
        pass
    else:
        raise ValueError(str(wildcards))
    return base_hrs * attempt


def set_winnowmap_preset(wildcards):
    preset = None
    if 'MAP-TO_HIFIEC' in wildcards.readset:
        preset = 'map-pb'
    elif 'ONTUL' in wildcards.readset:
        preset = 'map-ont'
    elif 'HIFIEC' in wildcards.readset:
        preset = 'map-pb'
    elif 'HIFIAF' in wildcards.readset:
        preset = 'map-pb'
    else:
        raise ValueError(str(wildcards))
    assert preset is not None
    return preset


rule qc_wmap_align_readsets:
    input:
        reads = select_winnowmap_reads,
        reference = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta'),
        ref_repkmer = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.k15.rep-grt09998.txt')
    output:
        bam = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.psort.bam',
        bai = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.psort.bam.bai'
    log:
        'log/output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.wmap.log'
    benchmark:
        'rsrc/output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.wmap.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: set_winnowmap_memory(wildcards, attempt),
        runtime_hrs = lambda wildcards, attempt: set_winnowmap_runtime(wildcards, attempt),
        mem_sort_mb = 4096,
        align_threads = config['num_cpu_high'] - config['num_cpu_low'],
        sort_threads = config['num_cpu_low'],
    params:
        individual = lambda wildcards: wildcards.sample,
        readgroup_id = lambda wildcards: wildcards.readset.replace('.', ''),
        preset = lambda wildcards: set_winnowmap_preset(wildcards),
        temp_prefix = lambda wildcards: f'temp/wmap/{wildcards.reference}/{wildcards.sample}/{wildcards.readset}/tmp_stsort_',
        temp_dir = lambda wildcards: f'temp/wmap/{wildcards.reference}/{wildcards.sample}/{wildcards.readset}/',
    shell:
        'rm -rfd {params.temp_dir} && mkdir -p {params.temp_dir} && '
        'winnowmap -W {input.ref_repkmer} -k 15 -t {resources.align_threads} -Y -L --eqx --MD -a -x {params.preset} '
        '-R "@RG\\tID:{params.readgroup_id}\\tSM:{params.individual}" --secondary=no '
        '{input.reference} {input.reads} | '
        'samtools view -u -F 260 | '
        'samtools sort -m {resources.mem_sort_mb}M -@ {resources.sort_threads} -T {params.temp_prefix} -O BAM > {output.bam} '
        ' && '
        'samtools index {output.bam} ; rm -rfd {params.temp_dir}'


rule dump_alignment_coverage:
    input:
        bam = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.psort.bam',
        bai = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.psort.bam.bai',
    output:
        bed = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.aln.bed',
    benchmark:
        'rsrc/output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.aln.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 1
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt,
    shell:
        'bedtools bamtobed -i {input.bam} > {output.bed}'
