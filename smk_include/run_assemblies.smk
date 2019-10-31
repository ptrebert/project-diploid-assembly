
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'

localrules: derive_wtdbg_parameter_preset, derive_flye_parameter_preset


rule derive_wtdbg_parameter_preset:
    input:
        '{filepath}.fastq.gz'
    output:
        '{filepath}.preset.wtdbg'
    run:
        import os
        preset = None
        file_name = os.path.basename(input[0])
        _, _, platform_spec = file_name.split('_')[:3]
        if platform_spec.startswith('pb'):
            if 'ccs' in platform_spec:
                preset = 'ccs'
            elif platform_spec.startswith('pbsq'):
                preset = 'sq'
            else:
                raise ValueError('Cannot handle PacBio platform spec: {} / {}'.format(file_name, platform_spec))
        elif platform_spec.startswith('ont'):
            preset = 'ont'
        else:
            raise ValueError('Unexpected platform spec: {} / {}'.format(file_name, platform_spec))

        if preset is None:
            raise ValueError('wtdbg preset is None: {}'.format(input[0]))

        with open(output[0], 'w') as dump:
            _ = dump.write(preset)


rule derive_flye_parameter_preset:
    input:
        '{filepath}.fastq.gz'
    output:
        '{filepath}.preset.flye'
    run:
        import os
        preset = None
        file_name = os.path.basename(input[0])
        _, _, platform_spec = file_name.split('_')[:3]
        if platform_spec.startswith('pb'):
            if 'ccs' in platform_spec:
                preset = '--pacbio-corr'
            elif 'clr' in platform_spec:
                preset = '--pacbio-raw'
            else:
                raise ValueError('Cannot handle PacBio platform spec: {} / {}'.format(file_name, platform_spec))
        elif platform_spec.startswith('ont'):
            preset = '--nano-raw'
        else:
            raise ValueError('Unexpected platform spec: {} / {}'.format(file_name, platform_spec))

        if preset is None:
            raise ValueError('flye preset is None: {}'.format(input[0]))

        with open(output[0], 'w') as dump:
            _ = dump.write(preset)


rule compute_wtdbg_squashed_assembly_layout:
    """
    Squashed assembly = replaces known reference (such as hg38)
    """
    input:
        fastq = 'input/fastq/complete/{sample}_1000.fastq.gz',
        preset = 'input/fastq/complete/{sample}_1000.preset.wtdbg'
    output:
        layout = 'references/assemblies/squashed_layout/wtdbg2/{sample}/{sample}.ctg.lay.gz',
    log: 'log/references/assemblies/squashed_layout/wtdbg2/{sample}.layout.log',
    benchmark: 'run/references/assemblies/squashed_layout/wtdbg2/{sample}.layout.rsrc',
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 409600
    params:
        param_preset = load_preset_file,
        out_prefix = lambda wildcards, output: output.layout.rsplit('.', 3)[0]
    shell:
        'wtdbg2 -x {params.param_preset} -i {input.fastq} -g3g -t {threads}' \
            ' -o {params.out_prefix} &> {log}'


rule compute_wtdbg_squashed_assembly_consensus:
    input:
        layout = 'references/assemblies/squashed_layout/wtdbg2/{sample}/{sample}.ctg.lay.gz'
    output:
        squashed_assembly = 'references/assemblies/{sample}_sqa-wtdbg.fasta'
    log: 'log/references/assemblies/{sample}_sqa-wtdbg.consensus.log'
    benchmark: 'run/references/assemblies/{sample}_sqa-wtdbg.consensus.rsrc'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 384,
        mem_total_mb = 12288
    shell:
        'wtpoa-cns -t {threads} -i {input.layout} -o {output.squashed_assembly} &> {log}'


rule compute_flye_squashed_assembly_layout:
    """
    Squashed assembly = replaces known reference (such as hg38)
    """
    input:
        fastq = 'input/fastq/complete/{sample}_1000.fastq.gz',
        preset = 'input/fastq/complete/{sample}_1000.preset.flye'
    output:
        layout = directory('references/assemblies/squashed_layout/flye/{sample}/00-assembly'),
        consensus = directory('references/assemblies/squashed_layout/flye/{sample}/10-consensus'),
        repeat = directory('references/assemblies/squashed_layout/flye/{sample}/20-repeat'),
        contigger = directory('references/assemblies/squashed_layout/flye/{sample}/30-contigger'),
        polishing = directory('references/assemblies/squashed_layout/flye/{sample}/40-polishing'),
        assm_graph = 'references/assemblies/squashed_layout/flye/{sample}/assembly_graph.gfa',
        assm_gv = 'references/assemblies/squashed_layout/flye/{sample}/assembly_graph.gv',
        assm_info = 'references/assemblies/squashed_layout/flye/{sample}/assembly_info.txt',
        run_params = 'references/assemblies/squashed_layout/flye/{sample}/params.json',
        assembly = 'references/assemblies/squashed_layout/flye/{sample}/assembly.fasta',
        final = 'references/assemblies/{sample}_sqa-flye.fasta',
    log: 'log/references/assemblies/squashed_layout/flye/{sample}.layout.log',
    benchmark: 'run/references/assemblies/squashed_layout/flye/{sample}.layout.rsrc',
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = int(1048576 / int(config['num_cpu_high'])),
        mem_total_mb = 1048576
    params:
        param_preset = load_preset_file,
        #param_preset = lambda wildcards, input: open(input.preset).read().strip(),
        out_prefix = lambda wildcards, output: os.path.dirname(output.assembly)
    shell:
        'flye {params.param_preset} {input.fastq} -g3g -t {threads}' \
            ' --debug --out-dir {params.out_prefix} &> {log}' \
            ' && ' \
            'cp {output.assembly} {output.final}'
