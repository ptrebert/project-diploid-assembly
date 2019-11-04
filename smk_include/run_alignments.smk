
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'canonical_dga.smk'
include: 'strandseq_dga_joint.smk'
include: 'strandseq_dga_split.smk'

localrules: master_run_alignments

rule master_run_alignments:
    input:


rule derive_minimap_parameter_preset:
    input:
        '{filepath}.fastq.gz'
    output:
        '{filepath}.preset.minimap'
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
        'run/output/alignments/reads_to_reference/{folder_path}/{sample}_map-to_{reference}.rsrc'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 2048,
        mem_total_mb = 65536
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
        'run/output/alignments/reads_to_reference/{folder_path}/{sample}_map-to_{reference}.pbn.rsrc'
    conda:
        config['conda_env_pbtools']
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 3072,
        mem_total_mb = 94208,
    params:
        align_threads = config['num_cpu_high'] - config['num_cpu_low'],
        sort_threads = config['num_cpu_low'],
        sort_memory_mb = 3072,  # also per thread
        preset = load_preset_file,
        individual = lambda wildcards: wildcards.sample.split('_')[0]
    shell:
        'pbmm2 align --log-level INFO --sort --sort-memory {params.sort_memory}M --no-bai ' \
            ' --alignment-threads {params.align_threads} --sort-threads {params.sort_threads} ' \
            ' --preset {params.preset} --sample {params.individual} ' \
            ' {input.reference} {input.reads} {output.bam} &> {log}'


rule bwa_strandseq_to_reference_alignment:
    input:
        mate1 = 'input/fastq/strand-seq/{individual}_{bioproject}/{individual}_{sample_id}_1.fastq.gz',
        mate2 = 'input/fastq/strand-seq/{individual}_{bioproject}/{individual}_{sample_id}_2.fastq.gz',
        ref_index = 'output/reference_assembly/squashed/bwa_index/{reference}.bwt',
    output:
        bam = temp('output/alignments/strandseq_to_reference/{reference}/{bioproject}/temp/aln/{individual}_{sample_id}.filt.sam.bam')
    log:
        bwa = 'log/output/alignments/strandseq_to_reference/{reference}/{bioproject}/{individual}_{sample_id}.bwa.log',
        samtools = 'log/output/alignments/strandseq_to_reference/{reference}/{bioproject}/{individual}_{sample_id}.samtools.log',
    benchmark:
        'run/output/alignments/strandseq_to_reference/{reference}/{bioproject}/{individual}_{sample_id}.rsrc'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = 2048,
        mem_total_mb = config['num_cpu_low'] * 2048
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
        reads = 'output/diploid_assembly/{folder_path}/draft/haploid_fastq/{subfolder}{readset}.fastq.gz',
        preset = 'output/diploid_assembly/{folder_path}/draft/haploid_fastq/{subfolder}{readset}.preset.minimap',
        contigs = 'output/diploid_assembly/{folder_path}/draft/haploid_fasta/{subfolder}{assembly}.fasta',
    output:
        sam = 'output/diploid_assembly/{folder_path}/polishing/alignments/{subfolder}{readset}_map-to_{assembly}.racon-p1.psort.sam',
    log:
        'log/output/diploid_assembly/{folder_path}/polishing/alignments/{subfolder}{readset}_map-to_{assembly}.racon-p1.psort.log',
    benchmark:
        'run/output/diploid_assembly/{folder_path}/polishing/alignments/{subfolder}{readset}_map-to_{assembly}.racon-p1.psort.rsrc',
    wildcard_constraints:
        subfolder = '(^$|splits\/)'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 768,
        mem_total_mb = 25600
    params:
        individual = lambda wildcards: wildcards.readset.split('_')[0],
        preset = load_preset_file
    shell:
        'minimap2 -t {threads} {params.preset} -R "@RG\\tID:1\\tSM:{params.individual}" ' \
            ' {input.contigs} {input.reads} 2> {log} | samtools sort | ' \
            ' samtools view -q 10 -F0x704 /dev/stdin > {output.sam}'


rule pbmm2_arrow_polish_alignment_pass1:
    input:
        reads = 'output/diploid_assembly/{folder_path}/draft/haploid_bam/{subfolder}{readset}.pbn.bam',
        preset = 'output/diploid_assembly/{folder_path}/draft/haploid_bam/{subfolder}{readset}.preset.pbmm2',
        contigs = 'output/diploid_assembly/{folder_path}/draft/haploid_fasta/{subfolder}{assembly}.fasta',
    output:
        bam = 'output/diploid_assembly/{folder_path}/polishing/alignments/{subfolder}{readset}_map-to_{assembly}.arrow-p1.psort.pbn.bam',
    log:
        'log/output/diploid_assembly/{folder_path}/polishing/alignments/{subfolder}{readset}_map-to_{assembly}.arrow-p1.psort.log',
    benchmark:
        'run/output/diploid_assembly/{folder_path}/polishing/alignments/{subfolder}{readset}_map-to_{assembly}.arrow-p1.psort.rsrc',
    wildcard_constraints:
        subfolder = '(^$|splits\/)'
    conda:
        config['conda_env_pbtools']
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 3072,
        mem_total_mb = 94208,
    params:
        align_threads = config['num_cpu_high'] - config['num_cpu_low'],
        sort_threads = config['num_cpu_low'],
        sort_memory_mb = 3072,  # also per thread
        preset = load_preset_file,
        individual = lambda wildcards: wildcards.readset.split('_')[0]
    shell:
        'pbmm2 align --log-level INFO --sort --sort-memory {params.sort_memory} --no-bai ' \
            ' --alignment-threads {params.align_threads} --sort-threads {params.sort_threads} ' \
            ' --preset {params.param_preset} --sample {params.individual} ' \
            ' {input.contigs} {input.reads} {output.bam} &> {log}'
