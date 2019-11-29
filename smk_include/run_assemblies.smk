
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'canonical_dga.smk'
include: 'strandseq_dga_joint.smk'
include: 'strandseq_dga_split.smk'


rule derive_wtdbg_parameter_preset:
    input:
        '{filepath}.fastq.gz'
    output:
        '{filepath}.preset.wtdbg'
    resources:
        runtime_hrs = 0,
        runtime_min = 10
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
                raise ValueError('Cannot handle PacBio platform spec: {threads} / {threads}'.format(file_name, platform_spec))
        elif platform_spec.startswith('ont'):
            preset = 'ont'
        else:
            raise ValueError('Unexpected platform spec: {threads} / {threads}'.format(file_name, platform_spec))

        if preset is None:
            raise ValueError('wtdbg preset is None: {threads}'.format(input[0]))

        with open(output[0], 'w') as dump:
            _ = dump.write(preset)


rule derive_flye_parameter_preset:
    input:
        '{filepath}.fastq.gz'
    output:
        '{filepath}.preset.flye'
    resources:
        runtime_hrs = 0,
        runtime_min = 10
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
                raise ValueError('Cannot handle PacBio platform spec: {threads} / {threads}'.format(file_name, platform_spec))
        elif platform_spec.startswith('ont'):
            preset = '--nano-raw'
        else:
            raise ValueError('Unexpected platform spec: {threads} / {threads}'.format(file_name, platform_spec))

        if preset is None:
            raise ValueError('flye preset is None: {threads}'.format(input[0]))

        with open(output[0], 'w') as dump:
            _ = dump.write(preset)


rule compute_wtdbg_nonhapres_assembly_layout:
    """
    Non-haplotype resolved assembly = replaces known reference (such as hg38)
    """
    input:
        fastq = 'input/fastq/complete/{sample}.fastq.gz',
        preset = 'input/fastq/complete/{sample}.preset.wtdbg'
    output:
        layout = 'output/reference_assembly/non-hap-res/layout/wtdbg2/{sample}/{sample}.ctg.lay.gz',
        aux = expand('output/reference_assembly/non-hap-res/layout/wtdbg2/{{sample}}/{{sample}}.{ext}',
                      ext=['1.dot.gz', '1.nodes', '1.reads', '2.dot.gz', '3.dot.gz',
                            'alignments.gz', 'ctg.dot.gz', 'events', 'frg.dot.gz', 'frg.nodes'])
    log:
        'log/output/reference_assembly/non-hap-res/layout/wtdbg2/{sample}_nhr-wtdbg.layout.log',
    benchmark:
        'run/output/reference_assembly/non-hap-res/layout/wtdbg2/{{sample}}_nhr-wtdbg.layout.t{}.rsrc'.format(config['num_cpu_max'])
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = lambda wildcards: int((163840 if '-ccs' in wildcards.sample else 393216) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards: 163840 if '-ccs' in wildcards.sample else 393216,
        runtime_hrs = lambda wildcards: 12 if '-ccs' in wildcards.sample else 36
    params:
        param_preset = load_preset_file,
        out_prefix = lambda wildcards, output: output.layout.rsplit('.', 3)[0]
    shell:
        'wtdbg2 -x {params.param_preset} -i {input.fastq} -g3g -t {threads}' \
            ' -o {params.out_prefix} &> {log}'


rule compute_wtdbg_nonhapres_assembly_consensus:
    input:
        layout = rules.compute_wtdbg_nonhapres_assembly_layout.output.layout
    output:
        nhr_assembly = protected('output/reference_assembly/non-hap-res/{sample}_nhr-wtdbg.fasta')
    log:
        'log/output/reference_assembly/non-hap-res/{sample}_nhr-wtdbg.consensus.log'
    benchmark:
        'run/output/reference_assembly/non-hap-res/{{sample}}_nhr-wtdbg.consensus.t{}.rsrc'.format(config['num_cpu_high'])
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards: int((8192 if '-ccs' in wildcards.sample else 24576) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards: 8192 if '-ccs' in wildcards.sample else 24576,
        runtime_hrs = lambda wildcards: 4 if '-ccs' in wildcards.sample else 12
    shell:
        'wtpoa-cns -t {threads} -i {input.layout} -o {output.nhr_assembly} &> {log}'


rule compute_flye_nonhapres_assembly:
    """
    Non-haplotype resolved assembly = replaces known reference (such as hg38)
    """
    input:
        fastq = 'input/fastq/complete/{sample}.fastq.gz',
        preset = 'input/fastq/complete/{sample}.preset.flye'
    output:
        layout = directory('output/reference_assembly/non-hap-res/layout/flye/{sample}/00-assembly'),
        consensus = directory('output/reference_assembly/non-hap-res/layout/flye/{sample}/10-consensus'),
        repeat = directory('output/reference_assembly/non-hap-res/layout/flye/{sample}/20-repeat'),
        contigger = directory('output/reference_assembly/non-hap-res/layout/flye/{sample}/30-contigger'),
        polishing = directory('output/reference_assembly/non-hap-res/layout/flye/{sample}/40-polishing'),
        assm_graph = 'output/reference_assembly/non-hap-res/layout/flye/{sample}/assembly_graph.gfa',
        assm_gv = 'output/reference_assembly/non-hap-res/layout/flye/{sample}/assembly_graph.gv',
        assm_info = 'output/reference_assembly/non-hap-res/layout/flye/{sample}/assembly_info.txt',
        run_params = 'output/reference_assembly/non-hap-res/layout/flye/{sample}/params.json',
        assm_source = 'output/reference_assembly/non-hap-res/layout/flye/{sample}/assembly.fasta',
        assembly = protected('output/reference_assembly/non-hap-res/{sample}_nhr-flye.fasta'),
    log:
        'log/output/reference_assembly/non-hap-res/{sample}_nhr-flye.layout.log',
    benchmark:
        'run/output/reference_assembly/non-hap-res/{{sample}}_nhr-flye.layout.t{}.rsrc'.format(config['num_cpu_max'])
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = lambda wildcards: int((557056 if '-ccs' in wildcards.sample else 1433600) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards: 557056 if '-ccs' in wildcards.sample else 1433600,
        runtime_hrs = lambda wildcards: 30 if '-ccs' in wildcards.sample else 128
    params:
        param_preset = load_preset_file,
        out_prefix = lambda wildcards, output: os.path.dirname(output.assm_source)
    shell:
        'flye {params.param_preset} {input.fastq} -g3g -t {threads}' \
            ' --debug --out-dir {params.out_prefix} &> {log}' \
            ' && ' \
            'cp {output.assm_source} {output.assembly}'


rule write_peregrine_nonhapres_fofn:
    input:
        'input/fastq/complete/{sample}.fastq.gz'
    output:
        'input/fastq/complete/{sample}.fofn'
    resources:
        runtime_hrs = 0,
        runtime_min = 5
    params:
        mount_point = '/wd'
    run:
        with open(output[0], 'w') as fofn:
            _ = fofn.write(params.mount_point + '/' + input[0])


rule compute_peregrine_nonhapres_assembly:
    input:
        container = 'output/container/docker/cschin/peregrine_{}.sif'.format(config['peregrine_version']),
        fofn = 'input/fastq/complete/{sample}.fofn'
    output:
        dir_seqdb = directory('output/reference_assembly/non-hap-res/layout/peregrine/{sample}/0-seqdb'),
        dir_index = directory('output/reference_assembly/non-hap-res/layout/peregrine/{sample}/1-index'),
        dir_ovlp = directory('output/reference_assembly/non-hap-res/layout/peregrine/{sample}/2-ovlp'),
        dir_asm = directory('output/reference_assembly/non-hap-res/layout/peregrine/{sample}/3-asm'),
        dir_cns = directory('output/reference_assembly/non-hap-res/layout/peregrine/{sample}/4-cns'),
        linked_asm = 'output/reference_assembly/non-hap-res/layout/peregrine/{sample}/p_ctg_cns.fa',
        assm = 'output/reference_assembly/non-hap-res/{sample}_nhr-pereg.fasta'
    log:
        pereg = 'log/output/reference_assembly/non-hap-res/{sample}_nhr-pereg.log',
        copy = 'log/output/reference_assembly/non-hap-res/{sample}_nhr-pereg.copy.log'
    benchmark:
        'run/output/reference_assembly/non-hap-res/{{sample}}_nhr-pereg.t{}.rsrc'.format(config['num_cpu_max'])
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = int(524288 / config['num_cpu_max']),
        mem_total_mb = 524288,
        runtime_hrs = 20
    params:
        bind_folder = lambda wildcards: os.getcwd(),
        out_folder = lambda wildcards, output: os.path.split(output.dir_seqdb)[0]
    shell:
        'yes yes | singularity run --bind {params.bind_folder}:/wd {input.container} ' \
            ' asm {input.fofn} {threads} {threads} {threads} {threads} {threads} {threads} {threads} {threads} {threads} ' \
            ' --with-consensus --shimmer-r 3 --best_n_ovlp 8 --output /wd/{params.out_folder} &> {log.pereg} ' \
            ' && '
            ' cp {output.linked_asm} {output.assm} &> {log.copy}'


### ----- Below this point ------ ###
### Haplotype-resolved assemblies ###


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
    log:
        'log/output/diploid_assembly/{variant}/{folder_path}/draft/temp/layout/wtdbg2/{hap_reads}.{hap}.wtdbg-layout.log'
    benchmark:
        'run/output/diploid_assembly/{{variant}}/{{folder_path}}/draft/temp/layout/wtdbg2/{{hap_reads}}.{{hap}}.wtdbg-layout.t{}.rsrc'.format(config['num_cpu_max'])
    wildcard_constraints:
        variant = '(canonical|strandseq_joint)'
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = lambda wildcards: int((163840 if '-ccs' in wildcards.hap_reads else 393216) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards: 163840 if '-ccs' in wildcards.hap_reads else 393216,
        runtime_hrs = lambda wildcards: 12 if '-ccs' in wildcards.hap_reads else 36
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
    log:
        'log/output/diploid_assembly/{variant}/{folder_path}/draft/haploid_fasta/{hap_reads}.{hap}.wtdbg-consensus.log'
    benchmark:
        'run/output/diploid_assembly/{{variant}}/{{folder_path}}/draft/haploid_fasta/{{hap_reads}}.{{hap}}.wtdbg-consensus.t{}.rsrc'.format(config['num_cpu_high'])
    wildcard_constraints:
            variant = '(canonical|strandseq_joint)'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards: int((8192 if '-ccs' in wildcards.hap_reads else 24576) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards: 8192 if '-ccs' in wildcards.hap_reads else 24576,
        runtime_hrs = lambda wildcards: 4 if '-ccs' in wildcards.hap_reads else 12
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
    log:
        'log/output/diploid_assembly/{variant}/{folder_path}/draft/haploid_fasta/{hap_reads}.{hap}.flye.log'
    benchmark:
        'run/output/diploid_assembly/{{variant}}/{{folder_path}}/draft/haploid_fasta/{{hap_reads}}.{{hap}}.flye.t{}.rsrc'.format(config['num_cpu_max'])
    wildcard_constraints:
            variant = '(canonical|strandseq_joint)'
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = lambda wildcards: int((557056 if '-ccs' in wildcards.hap_reads else 1433600) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards: 557056 if '-ccs' in wildcards.hap_reads else 1433600,
        runtime_hrs = lambda wildcards: 30 if '-ccs' in wildcards.hap_reads else 128
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
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/wtdbg2/{hap_reads}.{hap}.{sequence}.wtdbg-layout.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT_PROTECTED + '/draft/temp/layout/wtdbg2/{{hap_reads}}.{{hap}}.{{sequence}}.wtdbg-layout.t{}.rsrc'.format(config['num_cpu_high'])
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((4096 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: 1 * attempt
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
        haploid_assembly = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-wtdbg.{hap}.{sequence}.fasta'
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}.{hap}.{sequence}.wtdbg-consensus.log'
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT_PROTECTED + '/draft/haploid_fasta/{{hap_reads}}.{{hap}}.{{sequence}}.wtdbg-consensus.t{}.rsrc'.format(config['num_cpu_high'])
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((4096 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: 1 * attempt
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
        assembly = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-flye.{hap}.{sequence}.fasta',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}.{hap}.{sequence}.flye.log'
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT_PROTECTED + '/draft/haploid_fasta/{{hap_reads}}.{{hap}}.{{sequence}}.flye.t{}.rsrc'.format(config['num_cpu_high'])
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards: int((24576 if '-ccs' in wildcards.hap_reads else 131072) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards: 24576 if '-ccs' in wildcards.hap_reads else 131072,
        runtime_hrs = lambda wildcards: 0 if '-ccs' in wildcards.hap_reads else 5
    params:
        param_preset = load_preset_file,
        seq_len = load_seq_length_file,
        out_prefix = lambda wildcards, output: os.path.dirname(output.assm_source)
    shell:
        'flye {params.param_preset} {input.fastq} -g {params.seq_len} -t {threads}' \
            ' --debug --out-dir {params.out_prefix} &> {log}' \
            ' && ' \
            'cp {output.assm_source} {output.assembly}'
