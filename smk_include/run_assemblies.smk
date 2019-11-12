
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'canonical_dga.smk'
include: 'strandseq_dga_joint.smk'
include: 'strandseq_dga_split.smk'

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
        layout = 'output/reference_assembly/squashed/layout/wtdbg2/{sample}/{sample}.ctg.lay.gz',
        aux = expand('output/reference_assembly/squashed/layout/wtdbg2/{{sample}}/{{sample}}.{ext}',
                      ext=['1.dot.gz', '1.nodes', '1.reads', '2.dot.gz', '3.dot.gz',
                            'alignments.gz', 'binkmer', 'closed_bins', 'clps',
                            'ctg.dot.gz', 'events', 'frg.dot.gz', 'frg.nodes', 'kmerdep'])
    log: 'log/output/reference_assembly/squashed/layout/wtdbg2/{sample}.layout.log',
    benchmark: 'run/output/reference_assembly/squashed/layout/wtdbg2/{sample}.layout.rsrc',
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
        layout = rules.compute_wtdbg_squashed_assembly_layout.output.layout
    output:
        squashed_assembly = protected('output/reference_assembly/squashed/{sample}_sqa-wtdbg.fasta')
    log: 'log/output/reference_assembly/squashed/{sample}_sqa-wtdbg.consensus.log'
    benchmark: 'run/output/reference_assembly/squashed/{sample}_sqa-wtdbg.consensus.rsrc'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 384,
        mem_total_mb = 12288
    shell:
        'wtpoa-cns -t {threads} -i {input.layout} -o {output.squashed_assembly} &> {log}'


rule compute_flye_squashed_assembly:
    """
    Squashed assembly = replaces known reference (such as hg38)
    """
    input:
        fastq = 'input/fastq/complete/{sample}_1000.fastq.gz',
        preset = 'input/fastq/complete/{sample}_1000.preset.flye'
    output:
        layout = directory('output/reference_assembly/squashed/layout/flye/{sample}/00-assembly'),
        consensus = directory('output/reference_assembly/squashed/layout/flye/{sample}/10-consensus'),
        repeat = directory('output/reference_assembly/squashed/layout/flye/{sample}/20-repeat'),
        contigger = directory('output/reference_assembly/squashed/layout/flye/{sample}/30-contigger'),
        polishing = directory('output/reference_assembly/squashed/layout/flye/{sample}/40-polishing'),
        assm_graph = 'output/reference_assembly/squashed/layout/flye/{sample}/assembly_graph.gfa',
        assm_gv = 'output/reference_assembly/squashed/layout/flye/{sample}/assembly_graph.gv',
        assm_info = 'output/reference_assembly/squashed/layout/flye/{sample}/assembly_info.txt',
        run_params = 'output/reference_assembly/squashed/layout/flye/{sample}/params.json',
        assm_source = 'output/reference_assembly/squashed/layout/flye/{sample}/assembly.fasta',
        assembly = protected('output/reference_assembly/squashed/{sample}_sqa-flye.fasta'),
    log: 'log/output/reference_assembly/squashed/{sample}_sqa-flye.layout.log',
    benchmark: 'run/output/reference_assembly/squashed/{sample}_sqa-flye.layout.rsrc',
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = int(1048576 / int(config['num_cpu_high'])),
        mem_total_mb = 1048576
    params:
        param_preset = load_preset_file,
        out_prefix = lambda wildcards, output: os.path.dirname(output.assm_source)
    shell:
        'flye {params.param_preset} {input.fastq} -g3g -t {threads}' \
            ' --debug --out-dir {params.out_prefix} &> {log}' \
            ' && ' \
            'cp {output.assm_source} {output.assembly}'


rule compute_wtdbg_haploid_assembly_layout:
    input:
        fastq = 'output/diploid_assembly/{variant}/{folder_path}/draft/haploid_fastq/{hap_reads}.{hap}.fastq.gz',
        preset = 'output/diploid_assembly/{variant}/{folder_path}/draft/haploid_fastq/{hap_reads}.{hap}.preset.wtdbg'
    output:
        layout = 'output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/wtdbg2/{hap_reads}.{hap}/{hap_reads}.{hap}.ctg.lay.gz',
        # note: it seems that wtdbg does not always create all output files, so the list here
        # is incomplete - hopefully, it is enough to trigger a clean start if needed
        aux = expand('output/diploid_assembly/{{variant}}/{{folder_path}}/draft/temp/layout/wtdbg2/{{hap_reads}}.{{hap}}/{{hap_reads}}.{{hap}}.{ext}',
                      ext=['1.dot.gz', '1.nodes', '1.reads', '2.dot.gz', '3.dot.gz',
                            'alignments.gz', 'ctg.dot.gz', 'frg.dot.gz', 'frg.nodes'
                          ])
    log: 'log/output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/wtdbg2/{hap_reads}.{hap}.layout.log',
    benchmark: 'run/output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/wtdbg2/{hap_reads}.{hap}.layout.rsrc',
    wildcard_constraints:
        variant = '(canonical|strandseq_joint)'
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


rule compute_wtdbg_haploid_assembly_consensus:
    input:
        layout = rules.compute_wtdbg_haploid_assembly_layout.output.layout
    output:
        haploid_assembly = 'output/diploid_assembly/{variant}/{folder_path}/draft/haploid_fasta/{hap_reads}-wtdbg.{hap}.fasta'
    log: 'log/output/diploid_assembly/{variant}/{folder_path}/draft/haploid_fasta/{hap_reads}.{hap}.wtdbg.log'
    benchmark: 'run/output/diploid_assembly/{variant}/{folder_path}/draft/haploid_fasta/{hap_reads}.{hap}.wtdbg.rsrc'
    wildcard_constraints:
            variant = '(canonical|strandseq_joint)'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 384,
        mem_total_mb = 12288
    shell:
        'wtpoa-cns -t {threads} -i {input.layout} -o {output.haploid_assembly} &> {log}'


rule compute_flye_haploid_assembly:
    input:
        fastq = 'output/diploid_assembly/{variant}/{folder_path}/draft/haploid_fastq/{hap_reads}.{hap}.fastq.gz',
        preset = 'output/diploid_assembly/{variant}/{folder_path}/draft/haploid_fastq/{hap_reads}.{hap}.preset.flye'
    output:
        layout = directory('output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/flye/{hap_reads}.{hap}/00-assembly'),
        consensus = directory('output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/flye/{hap_reads}.{hap}/10-consensus'),
        repeat = directory('output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/flye/{hap_reads}.{hap}/20-repeat'),
        contigger = directory('output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/flye/{hap_reads}.{hap}/30-contigger'),
        polishing = directory('output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/flye/{hap_reads}.{hap}/40-polishing'),
        assm_graph = 'output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/flye/{hap_reads}.{hap}/assembly_graph.gfa',
        assm_gv = 'output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/flye/{hap_reads}.{hap}/assembly_graph.gv',
        assm_info = 'output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/flye/{hap_reads}.{hap}/assembly_info.txt',
        run_params = 'output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/flye/{hap_reads}.{hap}/params.json',
        assm_source = 'output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/flye/{hap_reads}.{hap}/assembly.fasta',
        assembly = 'output/diploid_assembly/{variant}/{folder_path}/draft/haploid_fasta/{hap_reads}-flye.{hap}.fasta',
    log: 'log/output/diploid_assembly/{variant}/{folder_path}/draft/haploid_fasta/{hap_reads}.{hap}.flye.log'
    benchmark: 'run/output/diploid_assembly/{variant}/{folder_path}/draft/haploid_fasta/{hap_reads}.{hap}.flye.rsrc'
    wildcard_constraints:
            variant = '(canonical|strandseq_joint)'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = int(1048576 / int(config['num_cpu_high'])),
        mem_total_mb = 1048576
    params:
        param_preset = load_preset_file,
        out_prefix = lambda wildcards, output: os.path.dirname(output.assm_source)
    shell:
        'flye {params.param_preset} {input.fastq} -g3g -t {threads}' \
            ' --debug --out-dir {params.out_prefix} &> {log}' \
            ' && ' \
            'cp {output.assm_source} {output.assembly}'

# Why do we need dedicated rules for computing the assemblies per haploid cluster/split?
# Setting the approx. assembly length requires accessing the respective reference
# FASTA index file

rule compute_wtdbg_haploid_split_assembly_layout:
    input:
        fastq = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq.gz',
        preset = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.preset.wtdbg',
        seq_info = 'output/reference_assembly/clustered/{sts_reads}/{reference}/sequences/{sequence}.seq'
    output:
        layout = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/wtdbg2/{hap_reads}.{hap}.{sequence}/{hap_reads}.{hap}.ctg.lay.gz',
        aux = expand('output/' + PATH_STRANDSEQ_DGA_SPLIT_PROTECTED + '/draft/temp/layout/wtdbg2/{{hap_reads}}.{{hap}}.{{sequence}}/{{hap_reads}}.{{hap}}.{ext}',
                      ext=['1.dot.gz', '1.nodes', '1.reads', '2.dot.gz', '3.dot.gz',
                            'alignments.gz', 'ctg.dot.gz', 'events', 'frg.dot.gz', 'frg.nodes'
                          ])
    log: 'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/wtdbg2/{hap_reads}.{hap}.{sequence}.layout.log',
    benchmark: 'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/wtdbg2/{hap_reads}.{hap}.{sequence}.layout.rsrc',
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 409600
    params:
        param_preset = load_preset_file,
        seq_len = load_seq_length_file,
        out_prefix = lambda wildcards, output: output.layout.rsplit('.', 3)[0]
    shell:
        'wtdbg2 -x {params.param_preset} -i {input.fastq} -g {params.seq_len} -t {threads}' \
            ' -o {params.out_prefix} &> {log}'


rule compute_wtdbg_haploid_split_assembly_consensus:
    input:
        layout = rules.compute_wtdbg_haploid_split_assembly_layout.output.layout
    output:
        haploid_assembly = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/split/{hap_reads}-wtdbg.{hap}.{sequence}.fasta'
    log: 'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}.{hap}.{sequence}.wtdbg.log'
    benchmark: 'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}.{hap}.{sequence}.wtdbg.rsrc'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 384,
        mem_total_mb = 12288
    shell:
        'wtpoa-cns -t {threads} -i {input.layout} -o {output.haploid_assembly} &> {log}'


rule compute_flye_haploid_split_assembly:
    input:
        fastq = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq.gz',
        preset = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.preset.flye',
        seq_info = 'output/reference_assembly/clustered/{sts_reads}/{reference}/sequences/{sequence}.seq'
    output:
        layout = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/flye/{hap_reads}.{hap}.{sequence}/00-assembly'),
        consensus = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/flye/{hap_reads}.{hap}.{sequence}/10-consensus'),
        repeat = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/flye/{hap_reads}.{hap}.{sequence}/20-repeat'),
        contigger = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/flye/{hap_reads}.{hap}.{sequence}/30-contigger'),
        polishing = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/flye/{hap_reads}.{hap}.{sequence}/40-polishing'),
        assm_graph = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/flye/{hap_reads}.{hap}.{sequence}/assembly_graph.gfa',
        assm_gv = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/flye/{hap_reads}.{hap}.{sequence}/assembly_graph.gv',
        assm_info = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/flye/{hap_reads}.{hap}.{sequence}/assembly_info.txt',
        run_params = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/flye/{hap_reads}.{hap}.{sequence}/params.json',
        assm_source = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/flye/{hap_reads}.{hap}.{sequence}/assembly.fasta',
        assembly = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/split/{hap_reads}-flye.{hap}.{sequence}.fasta',
    log: 'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}.{hap}.{sequence}.flye.log'
    benchmark: 'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}.{hap}.{sequence}.flye.rsrc'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = int(1048576 / int(config['num_cpu_high'])),
        mem_total_mb = 1048576
    params:
        param_preset = load_preset_file,
        seq_len = load_seq_length_file,
        out_prefix = lambda wildcards, output: os.path.dirname(output.assm_source)
    shell:
        'flye {params.param_preset} {input.fastq} -g {params.seq_len} -t {threads}' \
            ' --debug --out-dir {params.out_prefix} &> {log}' \
            ' && ' \
            'cp {output.assm_source} {output.assembly}'
