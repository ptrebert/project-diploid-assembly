
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'canonical_dga.smk'

localrules: master_run_alignments

rule master_run_alignments:
    input:


rule minimap_reads_to_reference_alignment:
    input:
        reads = 'input/fastq/complete/{sample}.fastq.gz',
        reference = 'references/assemblies/{reference}.fasta'
    output:
        'output/alignments/reads_to_reference/{sample}_map-to_{reference}.sam'
    log:
        'log/output/alignments/reads_to_reference/{sample}_map-to_{reference}.log'
    benchmark:
        'run/output/alignments/reads_to_reference/{sample}_map-to_{reference}.rsrc'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 2048,
        mem_total_mb = 65536
    run:
        preset = None
        is_pacbio = False
        individual = wildcards.sample.split('_')[0]
        if '_pb' in wildcards.sample:
            is_pacbio = True
            if '-clr' in wildcards.sample:
                preset = 'map-pb'
            elif '-ccs' in wildcards.sample:
                preset = 'asm20'
            else:
                raise ValueError('Unknown minimap2 preset: {}'.format(wildcards.sample))
        if '_ont' in wildcards.sample:
            preset = 'map-ont'

        assert preset is not None, 'No minimap preset selected for sample {}'.format(wildcards.sample)

        exec = 'minimap2 -t {threads} -x ' + preset
        if is_pacbio:
            exec += ' -H '
        exec += ' -a -o {output}'
        exec += ' -R "@RG\\tID:1\\tSM:' + individual + '" '
        exec += ' {input.reference} {input.reads}'
        exec += ' &> {log}'
        shell(exec)


rule pbmm2_reads_to_reference_alignment:
    input:
        reads = 'input/bam/complete/{sample}.pbn.bam',
        reference = 'references/assemblies/{reference}.fasta'
    output:
        bam = 'output/alignments/reads_to_reference/{sample}_map-to_{reference}.psort.pbn.bam',
    log:
        'log/output/alignments/reads_to_reference/{sample}_map-to_{reference}.pbn.log'
    benchmark:
        'run/output/alignments/reads_to_reference/{sample}_map-to_{reference}.pbn.rsrc'
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
        param_preset = lambda wildcards: config['pbmm2_presets'][wildcards.sample.rsplit('_', 1)[0]],
        individual = lambda wildcards: wildcards.sample.split('_')[0]
    shell:
        'pbmm2 align --log-level INFO --sort --sort-memory {params.sort_memory}M --no-bai ' \
            ' --alignment-threads {params.align_threads} --sort-threads {params.sort_threads} ' \
            ' --preset {params.param_preset} --sample {params.individual} ' \
            ' {input.reference} {input.reads} {output.bam} &> {log}'


rule bwa_strandseq_to_reference_alignment:
    input:
        mate1 = 'input/fastq/strand-seq/{individual}_{bioproject}/{individual}_{project}_{spec}_{lib_id}_{run_id}_1.fastq.gz',
        mate2 = 'input/fastq/strand-seq/{individual}_{bioproject}/{individual}_{project}_{spec}_{lib_id}_{run_id}_2.fastq.gz',
        ref_index = 'references/assemblies/bwa_index/{reference}.bwt',
    output:
        bam = 'output/alignments/strandseq_to_reference/{reference}.{individual}.{bioproject}/{individual}_{project}_{spec}_{lib_id}_{run_id}.filt.sam.bam'
    log:
        bwa = 'log/output/alignments/strandseq_to_reference/{reference}.{individual}.{bioproject}/{individual}_{project}_{spec}_{lib_id}_{run_id}.bwa.log',
        samtools = 'log/output/alignments/strandseq_to_reference/{reference}.{individual}.{bioproject}/{individual}_{project}_{spec}_{lib_id}_{run_id}.samtools.log',
    benchmark:
        'run/output/alignments/strandseq_to_reference/{reference}_{individual}_{bioproject}/{individual}_{project}_{spec}_{lib_id}_{run_id}.rsrc'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = 2048,
        mem_total_mb = config['num_cpu_low'] * 2048
    params:
        read_spec = lambda wildcards: wildcards.spec.split('-')[1],
        idx_prefix = lambda wildcards, input: input.ref_index.split('.')[0]
    shell:
        'bwa mem -t {threads}' \
            ' -R "@RG\\tID:{wildcards.individual}_{params.read_spec}\\tPL:Illumina\\tSM:{wildcards.individual}"' \
            ' -v 2 {params.idx_prefix} {input.mate1} {input.mate2} 2> {log.bwa} | ' \
            ' samtools view -b -F 2304 /dev/stdin > {output.bam} 2> {log.samtools}'


rule racon_strandseq_polish_alignment_pass1:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    pol_reads = FASTQ file used for Racon contig polishing
    """
    input:
        reads = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}-{assembler}/{vc_reads}/{sts_reads}/{pol_reads}.{hap}.fastq.gz',
        contigs = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}-{assembler}/{vc_reads}/{sts_reads}/consensus/{hap_reads}-{assembler}.{hap}.fasta'
    output:
        sam = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}-{assembler}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{reference}-{assembler}.{hap}.racon-p1.psort.sam',
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}-{assembler}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{reference}-{assembler}.{hap}.racon-p1.log'
    benchmark:
        'run/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}-{assembler}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{reference}-{assembler}.{hap}.racon-p1.rsrc'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 768,
        mem_total_mb = 25600
    run:
        preset = None
        is_pacbio = False
        individual = wildcards.pol_reads.split('_')[0]
        if '_pb' in wildcards.pol_reads:
            is_pacbio = True
            if '-clr' in wildcards.pol_reads:
                preset = 'map-pb'
            elif '-ccs' in wildcards.pol_reads:
                preset = 'asm20'
            else:
                raise ValueError('Unknown minimap2 preset: {}'.format(wildcards.pol_reads))
        if '_ont' in wildcards.pol_reads:
            preset = 'map-ont'

        assert preset is not None, 'No minimap preset selected for sample {}'.format(wildcards.pol_reads)

        exec = 'minimap2 -t {threads} -x ' + preset
        if is_pacbio:
            exec += ' -H '
        exec += ' -a'
        exec += ' --eqx -m 5000 --secondary=no'  # parameters recommended by Aaron Wenger
        exec += ' -R "@RG\\tID:1\\tSM:' + individual + '" '
        exec += ' {input.contigs} {input.reads}'
        exec += ' 2> {log}'
        exec += ' | samtools sort'
        exec += ' | samtools view -q 10 -F0x704 /dev/stdin'
        exec += ' > {output.sam}'
        shell(exec)


rule pbmm2_strandseq_polish_alignment_pass1:
    input:
        reads = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{pol_reads}.{hap}.pbn.bam',
        contigs = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}-{assembler}.{hap}.fasta'
    output:
        bam = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{sample}-{assembler}.{hap}.arrow-p1.psort.pbn.bam',
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{sample}-{assembler}.{hap}.arrow-p1.log',
    benchmark:
        'run/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{sample}-{assembler}.{hap}.arrow-p1.rsrc'
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
        param_preset = lambda wildcards: config['pbmm2_presets'][wildcards.pol_reads.rsplit('_', 1)[0]],
        individual = lambda wildcards: wildcards.pol_reads.split('_')[0]
    shell:
        'pbmm2 align --log-level INFO --sort --sort-memory {params.sort_memory} --no-bai ' \
            ' --alignment-threads {params.align_threads} --sort-threads {params.sort_threads} ' \
            ' --preset {params.param_preset} --sample {params.individual} ' \
            ' {input.contigs} {input.reads} {output.bam} &> {log}'
