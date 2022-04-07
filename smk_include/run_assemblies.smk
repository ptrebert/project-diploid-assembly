
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
        '{folder_path}/{file_name}.fastq.gz',
        '{folder_path}/{file_name}.stats'
    output:
        '{folder_path}/{file_name}.preset.flye'
    run:
        import os
        preset = None

        target_coverage = int(config['flye_target_coverage'])

        with open(input[1], 'r') as table:
            for line in table:
                if not line.startswith('cov_geq'):
                    continue
                rlen_info, cov = line.strip().split()  # this is the entry cov_geq_0
                cov = int(float(cov))
                if cov > target_coverage:
                    # only used for initial disjointig assembly
                    # for final assembly, all reads are used
                    preset = '--asm-coverage {} '.format(target_coverage)
                else:
                    preset = ''
                break

        file_name = os.path.basename(input[0])
        _, _, platform_spec = file_name.split('_')[:3]
        if platform_spec.startswith('pb'):
            if 'ccs' in platform_spec:
                preset += '--pacbio-hifi'
            elif 'clr' in platform_spec:
                preset += '--pacbio-raw'
            else:
                raise ValueError('Cannot handle PacBio platform spec: {} / {}'.format(file_name, platform_spec))
        elif platform_spec.startswith('ont'):
            preset += '--nano-raw'
        else:
            raise ValueError('Unexpected platform spec: {} / {}'.format(file_name, platform_spec))

        if preset is None:
            raise ValueError('flye preset is None: {}'.format(input[0]))

        with open(output[0], 'w') as dump:
            _ = dump.write(preset)


rule derive_canu_parameter_preset:
    input:
        '{folder_path}/{file_name}.fastq.gz',
        #'{folder_path}/{file_name}.stats'
    output:
        '{folder_path}/{file_name}.preset.canu'
    run:
        import os
        preset = ''

        file_name = os.path.basename(input[0])
        _, _, platform_spec = file_name.split('_')[:3]
        if platform_spec.startswith('pb'):
            if 'ccs' in platform_spec:
                preset += '-pacbio-hifi'
            elif 'clr' in platform_spec:
                preset += '-pacbio-raw'
            else:
                raise ValueError('Cannot handle PacBio platform spec: {} / {}'.format(file_name, platform_spec))
        elif platform_spec.startswith('ont'):
            preset += '-nano-raw'
        else:
            raise ValueError('Unexpected platform spec: {} / {}'.format(file_name, platform_spec))

        if not preset:
            raise ValueError('canu preset is None: {}'.format(input[0]))

        with open(output[0], 'w') as dump:
            _ = dump.write(preset)


rule compute_wtdbg_nonhapres_assembly_layout:
    """
    Non-haplotype resolved assembly = replaces known reference (such as hg38)
    """
    input:
        fastq = 'input/fastq/{sample}.fastq.gz',
        preset = 'input/fastq/{sample}.preset.wtdbg',
        seq_info = 'references/assemblies/' + config['use_genome_size'] +'.fasta.fai'
    output:
        layout = 'output/reference_assembly/non-hap-res/layout/wtdbg2/{sample}/{sample}.ctg.lay.gz',
        aux = expand('output/reference_assembly/non-hap-res/layout/wtdbg2/{{sample}}/{{sample}}.{ext}',
                      ext=['1.dot.gz', '1.nodes', '1.reads', '2.dot.gz', '3.dot.gz',
                            'alignments.gz', 'ctg.dot.gz', 'events', 'frg.dot.gz', 'frg.nodes'])
    log:
        'log/output/reference_assembly/non-hap-res/layout/wtdbg2/{sample}_nhr-wtdbg.layout.log',
    benchmark:
        'rsrc/output/reference_assembly/non-hap-res/layout/wtdbg2/{{sample}}_nhr-wtdbg.layout.t{}.rsrc'.format(config['num_cpu_max'])
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = lambda wildcards: int((163840 if '-ccs' in wildcards.sample else 393216) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards: 163840 if '-ccs' in wildcards.sample else 393216,
        runtime_hrs = lambda wildcards: 16 if '-ccs' in wildcards.sample else 40
    params:
        param_preset = load_preset_file,
        out_prefix = lambda wildcards, output: output.layout.rsplit('.', 3)[0],
        genome_size = load_seq_length_file
    shell:
        'wtdbg2 -x {params.param_preset} -i {input.fastq} -g {params.genome_size} -t {threads}'
            ' -o {params.out_prefix} &> {log}'


rule compute_wtdbg_nonhapres_assembly_consensus:
    input:
        layout = rules.compute_wtdbg_nonhapres_assembly_layout.output.layout
    output:
        nhr_assembly = 'output/reference_assembly/non-hap-res/{sample}_nhr-wtdbg.fasta'
    log:
        'log/output/reference_assembly/non-hap-res/{sample}_nhr-wtdbg.consensus.log'
    benchmark:
        'rsrc/output/reference_assembly/non-hap-res/{{sample}}_nhr-wtdbg.consensus.t{}.rsrc'.format(config['num_cpu_high'])
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards: int((8192 if '-ccs' in wildcards.sample else 24576) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards: 8192 if '-ccs' in wildcards.sample else 24576,
        runtime_hrs = lambda wildcards: 8 if '-ccs' in wildcards.sample else 20
    shell:
        'wtpoa-cns -t {threads} -i {input.layout} -o {output.nhr_assembly} &> {log}'


rule compute_flye_nonhapres_assembly:
    """
    Non-haplotype resolved assembly = replaces known reference (such as hg38)
    """
    input:
        fastq = 'input/fastq/{sample}.fastq.gz',
        preset = 'input/fastq/{sample}.preset.flye',
        seq_info = 'references/assemblies/' + config['use_genome_size'] +'.fasta.fai'
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
        assembly = 'output/reference_assembly/non-hap-res/{sample}_nhr-flye.fasta'
    log:
        'log/output/reference_assembly/non-hap-res/{sample}_nhr-flye.layout.log',
    benchmark:
        'rsrc/output/reference_assembly/non-hap-res/{{sample}}_nhr-flye.layout.t{}.rsrc'.format(config['num_cpu_max'])
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((1499136 if attempt <= 1 else 2949120) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards, attempt: 1499136 if attempt <= 1 else 2949120,
        runtime_hrs = lambda wildcards, attempt: 47 if '-ccs' in wildcards.sample else 119
    params:
        param_preset = load_preset_file,
        out_prefix = lambda wildcards, output: os.path.dirname(output.assm_source),
        genome_size = load_seq_length_file
    shell:
        'flye {params.param_preset} {input.fastq} -g {params.genome_size} -t {threads}'
            ' --debug --out-dir {params.out_prefix} &> {log}'
            ' && '
            'cp {output.assm_source} {output.assembly}'


def set_hifiasm_nhr_resources(wildcards):
    """
    regular HGSVC HiFi sample: ~30x
    regular HPRC HiFi sample: ~35x - ~40x
    high-coverage HGSVC sample: ~60x
    """
    if 'hprc' in wildcards.sample:
        threads = config['num_cpu_max'] // 2
        memory_mb = 229376
        runtime_hrs = 47
    elif wildcards.sample.startswith('HC'):
        # custom prefix for high-cov samples
        threads = config['num_cpu_max']
        memory_mb = 393216
        runtime_hrs = 59
    else:
        threads = config['num_cpu_high']
        memory_mb = 180224
        runtime_hrs = 47
    return threads, memory_mb, runtime_hrs


rule compute_hifiasm_nonhapres_assembly:
    """
    Runtime for slow I/O systems
    Memory consumption is fairly low and stable
    """
    input:
        fastq = 'input/fastq/{sample}.fastq.gz',
    output:
        primary_unitigs = 'output/reference_assembly/non-hap-res/layout/hifiasm/{sample}/{sample}.p_utg.gfa',
        primary_contigs = 'output/reference_assembly/non-hap-res/layout/hifiasm/{sample}/{sample}.p_ctg.gfa',
        raw_unitigs = 'output/reference_assembly/non-hap-res/layout/hifiasm/{sample}/{sample}.r_utg.gfa',
        discard = multiext(
                'output/reference_assembly/non-hap-res/layout/hifiasm/{sample}/{sample}',
                '.a_ctg.gfa', '.a_ctg.noseq.gfa',
                '.p_ctg.noseq.gfa', '.p_utg.noseq.gfa',
                '.r_utg.noseq.gfa'
            ),
        ec_reads = temp('output/reference_assembly/non-hap-res/layout/hifiasm/{sample}/{sample}.ec.fa'),
        reads_ava = temp('output/reference_assembly/non-hap-res/layout/hifiasm/{sample}/{sample}.ovlp.paf')
        # exclude the all-vs-all overlap information from output
        # => won't be deleted by snakemake in case of errors;
        # should save time if assembly job is restarted
        # '.ec.bin', '.ovlp.reverse.bin', '.ovlp.source.bin',
    log:
        hifiasm = 'log/output/reference_assembly/non-hap-res/{sample}_nhr-hifiasm.log',
    benchmark:
        os.path.join('rsrc/output/reference_assembly/non-hap-res',
                     '{sample}_nhr-hifiasm' + '.t{}.rsrc'.format(config['num_cpu_high']))
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: lambda wildcards: set_hifiasm_nhr_resources(wildcards)[0]
    resources:
        mem_total_mb = lambda wildcards, attempt: set_hifiasm_nhr_resources(wildcards)[1] * attempt,
        runtime_hrs = lambda wildcards, attempt: set_hifiasm_nhr_resources(wildcards)[2] * attempt
    params:
        prefix = lambda wildcards, output: output.primary_contigs.rsplit('.', 2)[0],
    shell:
        'hifiasm -o {params.prefix} -t {threads} --write-ec --write-paf --primary {input.fastq} &> {log.hifiasm}'


rule compress_hifiasm_ec_reads:
    input:
        'output/reference_assembly/non-hap-res/layout/hifiasm/{sample}/{sample}.ec.fa',
    output:
        'output/reference_assembly/non-hap-res/layout/hifiasm/{sample}.ec-reads.fasta.gz',
    conda:
        '../environment/conda/conda_shelltools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(1024 * attempt / config['num_cpu_low']),
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: 7 * attempt
    shell:
        'pigz -c --best -p {threads} {input} > {output}'


rule compress_hifiasm_reads_overlaps:
    input:
        'output/reference_assembly/non-hap-res/layout/hifiasm/{sample}/{sample}.ovlp.paf',
    output:
        'output/reference_assembly/non-hap-res/layout/hifiasm/{sample}.ava-reads.paf.gz',
    conda:
        '../environment/conda/conda_shelltools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(1024 * attempt / config['num_cpu_low']),
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'pigz -c --best -p {threads} {input} > {output}'


rule convert_nonhapres_gfa_to_fasta:
    """
    input includes gzipped error-corrected reads and read self-overlaps just to force creation
    """
    input:
        gfa = 'output/reference_assembly/non-hap-res/layout/hifiasm/{sample}/{sample}.p_ctg.gfa',
        paf_gz = 'output/reference_assembly/non-hap-res/layout/hifiasm/{sample}.ava-reads.paf.gz',
        fasta_gz = 'output/reference_assembly/non-hap-res/layout/hifiasm/{sample}.ec-reads.fasta.gz',
    output:
        fasta = protected('output/reference_assembly/non-hap-res/{sample}_nhr-hifiasm.fasta'),
        rc_map = 'output/reference_assembly/non-hap-res/{sample}_nhr-hifiasm.read-contig.map',
        stats = 'output/reference_assembly/non-hap-res/{sample}_nhr-hifiasm.contig.stats',
    log:
        'log/output/reference_assembly/non-hap-res/{sample}_nhr-hifiasm.gfa-convert.log'
    benchmark:
        'rsrc/output/reference_assembly/non-hap-res/{sample}_nhr-hifiasm.gfa-convert' + '.t{}.rsrc'.format(config['num_cpu_low'])
    conda:
        '../environment/conda/conda_pyscript.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = lambda wildcards: find_script_path('gfa_to_fasta.py')
    shell:
        '{params.script_exec} --gfa {input[0]} --n-cpus {threads} '
        '--out-fasta {output.fasta} --out-map {output.rc_map} --out-stats {output.stats}'


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
        'rsrc/output/diploid_assembly/{{variant}}/{{folder_path}}/draft/temp/layout/wtdbg2/{{hap_reads}}.{{hap}}.wtdbg-layout.t{}.rsrc'.format(config['num_cpu_max'])
    conda:
        '../environment/conda/conda_biotools.yml'
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
        'wtdbg2 -x {params.param_preset} -i {input.fastq} -g3g -t {threads}'
            ' -o {params.out_prefix} &> {log}'


rule compute_wtdbg_haploid_assembly_consensus:
    input:
        layout = rules.compute_wtdbg_haploid_assembly_layout.output.layout
    output:
        haploid_assembly = 'output/diploid_assembly/{variant}/{folder_path}/draft/haploid_assembly/{hap_reads}-wtdbg.{hap}.fasta'
    log:
        'log/output/diploid_assembly/{variant}/{folder_path}/draft/haploid_assembly/{hap_reads}.{hap}.wtdbg-consensus.log'
    benchmark:
        'rsrc/output/diploid_assembly/{{variant}}/{{folder_path}}/draft/haploid_assembly/{{hap_reads}}.{{hap}}.wtdbg-consensus.t{}.rsrc'.format(config['num_cpu_high'])
    conda:
        '../environment/conda/conda_biotools.yml'
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
        assembly = 'output/diploid_assembly/{variant}/{folder_path}/draft/haploid_assembly/{hap_reads}-flye.{hap}.fasta'
    log:
        'log/output/diploid_assembly/{variant}/{folder_path}/draft/haploid_assembly/{hap_reads}.{hap}.flye.log'
    benchmark:
        'rsrc/output/diploid_assembly/{{variant}}/{{folder_path}}/draft/haploid_assembly/{{hap_reads}}.{{hap}}.flye.t{}.rsrc'.format(config['num_cpu_max'])
    conda:
        '../environment/conda/conda_biotools.yml'
    wildcard_constraints:
            variant = '(canonical|strandseq_joint)'
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = lambda wildcards: int((557056 if '-ccs' in wildcards.hap_reads else 1433600) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards: 557056 if '-ccs' in wildcards.hap_reads else 1433600,
        runtime_hrs = lambda wildcards: 30 if '-ccs' in wildcards.hap_reads else 119
    params:
        param_preset = load_preset_file,
        out_prefix = lambda wildcards, output: os.path.dirname(output.assm_source)
    shell:
        'flye {params.param_preset} {input.fastq} -g3g -t {threads}'
            ' --debug --out-dir {params.out_prefix} &> {log}'
            ' && '
            'cp {output.assm_source} {output.assembly}'


# Why do we need dedicated rules for computing the assemblies per haploid cluster/split?
# Setting the approx. assembly length requires accessing the respective reference
# FASTA index file

rule compute_wtdbg_haploid_split_assembly_layout:
    input:
        fastq = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq.gz',
        preset = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.preset.wtdbg',
        seq_info = 'output/reference_assembly/clustered/{sseq_reads}/{reference}/sequences/{sequence}.seq'
    output:
        layout = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/wtdbg2/{hap_reads}.{hap}.{sequence}/{hap_reads}.{hap}.ctg.lay.gz',
        aux = expand('output/' + PATH_STRANDSEQ_DGA_SPLIT_PROTECTED + '/draft/temp/layout/wtdbg2/{{hap_reads}}.{{hap}}.{{sequence}}/{{hap_reads}}.{{hap}}.{ext}',
                      ext=['1.dot.gz', '1.nodes', '1.reads', '2.dot.gz', '3.dot.gz',
                            'alignments.gz', 'ctg.dot.gz', 'events', 'frg.dot.gz', 'frg.nodes'
                          ])
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/wtdbg2/{hap_reads}.{hap}.{sequence}.wtdbg-layout.log',
    benchmark:
        'rsrc/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/wtdbg2/{{hap_reads}}.{{hap}}.{{sequence}}.wtdbg-layout.t{}.rsrc'.format(config['num_cpu_high'])
    conda:
        '../environment/conda/conda_biotools.yml'
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
        'wtdbg2 -x {params.param_preset} -i {input.fastq} -g {params.seq_len} -t {threads}'
            ' -o {params.out_prefix} &> {log}'


rule compute_wtdbg_haploid_split_assembly_consensus:
    input:
        layout = rules.compute_wtdbg_haploid_split_assembly_layout.output.layout
    output:
        haploid_assembly = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-wtdbg.{hap}.{sequence}.fasta'
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}.{hap}.{sequence}.wtdbg-consensus.log'
    benchmark:
        'rsrc/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{{hap_reads}}.{{hap}}.{{sequence}}.wtdbg-consensus.t{}.rsrc'.format(config['num_cpu_high'])
    conda:
        '../environment/conda/conda_biotools.yml'
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
        seq_info = 'output/reference_assembly/clustered/{sseq_reads}/{reference}/sequences/{sequence}.seq'
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
        assembly = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-flye.{hap}.{sequence}.fasta',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}.{hap}.{sequence}.flye.log'
    benchmark:
        'rsrc/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{{hap_reads}}.{{hap}}.{{sequence}}.flye.t{}.rsrc'.format(config['num_cpu_high'])
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((110592 if attempt <= 1 else 188416 * (attempt - 1)) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 110592 if attempt <= 1 else 188416 * (attempt - 1),
        runtime_hrs = lambda wildcards, attempt: 8 * attempt
    params:
        param_preset = load_preset_file,
        seq_len = load_seq_length_file,
        out_prefix = lambda wildcards, output: os.path.dirname(output.assm_source)
    shell:
        'flye {params.param_preset} {input.fastq} -g {params.seq_len} -t {threads}'
            ' --debug --out-dir {params.out_prefix} &> {log}'
            ' && '
            'cp {output.assm_source} {output.assembly}'


rule compute_canu_haploid_split_assembly:
    """
    Because of the substantial resource requirements for Canu,
    this is currently only implemented for haploid split assemblies.
    """
    input:
        fastq = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq.gz',
        preset = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.preset.canu',
        seq_info = 'output/reference_assembly/clustered/{sseq_reads}/{reference}/sequences/{sequence}.seq'
    output:
        logs = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/canu/{hap_reads}.{hap}.{sequence}/canu-logs'),
        scripts = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/canu/{hap_reads}.{hap}.{sequence}/canu-scripts'),
        correct = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/canu/{hap_reads}.{hap}.{sequence}/correction'),
        trimming = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/canu/{hap_reads}.{hap}.{sequence}/trimming'),
        unitigging = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/canu/{hap_reads}.{hap}.{sequence}/unitigging'),
        # below: stuff which actually has the prefix
        seqstore = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/canu/{hap_reads}.{hap}.{sequence}/{hap_reads}.{hap}.{sequence}.seqStore'),
        stuff = multiext('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/canu/{hap_reads}.{hap}.{sequence}/{hap_reads}.{hap}.{sequence}',
                         '.contigs.layout', '.contigs.layout.readToTig', '.contigs.layout.tigInfo', '.unitigs.bed', '.unitigs.fasta', '.unitigs.gfa',
                         '.unitigs.layout', '.unitigs.layout.readToTig', '.unitigs.layout.tigInfo'),
        assm_source = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/canu/{hap_reads}.{hap}.{sequence}/{hap_reads}.{hap}.{sequence}.contigs.fasta',
        assembly = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-canu.{hap}.{sequence}.fasta'
    log:
        os.path.join('log/output', PATH_STRANDSEQ_DGA_SPLIT, 'draft/haploid_assembly/{hap_reads}.{hap}.{sequence}.canu.log')
    benchmark:
        os.path.join('rsrc/output', PATH_STRANDSEQ_DGA_SPLIT, 'draft/haploid_assembly/',
                     '{hap_reads}.{hap}.{sequence}.canu' + '.t{}.rsrc'.format(config['num_cpu_high']))
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 110592 if attempt <= 1 else 188416 * (attempt - 1),
        mem_per_cpu_mb = lambda wildcards, attempt: int((110592 if attempt <= 1 else 188416 * (attempt - 1)) / config['num_cpu_high']),
        runtime_hrs = lambda wildcards, attempt: 24 + 12 * attempt
    params:
        out_prefix = lambda wildcards: '.'.join([wildcards.hap_reads, wildcards.hap, wildcards.sequence]),
        out_dir = lambda wildcards, output: os.path.dirname(output.logs),
        param_preset = load_preset_file,
        seq_len = load_seq_length_file,
    shell:
        'canu -p {params.out_prefix} -d {params.out_dir} '
            'genomeSize={params.seq_len} '
            'useGrid=false saveOverlaps=false saveReadCorrections=false saveReads=false '
            '{params.param_preset} {input.fastq} &> {log} '
            ' && '
            'cp {output.assm_source} {output.assembly}'


rule reduce_haplotags_to_readlist:
    input:
        tags = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{hap_reads}.{sequence}.tags.fq.tsv',
    output:
        list_h1 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{hap_reads}.h1.{sequence}.tags.list',
        list_h2 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{hap_reads}.h2.{sequence}.tags.list'
    shell:
        'cut -f 1,2 {input} | egrep "H1" | cut -f 1 > {output.list_h1}'
        ' && '
        'cut -f 1,2 {input} | egrep "H2" | cut -f 1 > {output.list_h2}'


rule compute_hifiasm_haploid_split_assembly:
    """
    Runtime for slow I/O systems
    Memory consumption is fairly low and stable
    """
    input:
        fastq_h1 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.h1.{sequence}.fastq.gz',
        fastq_h2 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.h2.{sequence}.fastq.gz',
        fastq_un = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.un.{sequence}.fastq.gz',
        tags_h1 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{hap_reads}.h1.{sequence}.tags.list',
        tags_h2 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{hap_reads}.h2.{sequence}.tags.list',
    output:
        hap1_contigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}/{hap_reads}.{sequence}.dip.hap1.p_ctg.gfa',
        hap2_contigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}/{hap_reads}.{sequence}.dip.hap2.p_ctg.gfa',
        raw_unitigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}/{hap_reads}.{sequence}.dip.r_utg.gfa',
        discard = multiext(
                'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}/{hap_reads}.{sequence}.dip',
                '.hap1.p_ctg.noseq.gfa', '.hap2.p_ctg.noseq.gfa',
                '.hap1.p_ctg.lowQ.bed', '.hap2.p_ctg.lowQ.bed',
                '.p_utg.gfa', '.p_utg.lowQ.bed', '.p_utg.noseq.gfa', 
                '.r_utg.noseq.gfa', '.r_utg.lowQ.bed'
            )
        # as above: exclude read self overlaps from output
        # '.ec.bin', '.ovlp.reverse.bin', '.ovlp.source.bin',
    log:
        hifiasm = 'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}.hifiasm.log',
    benchmark:
        os.path.join('rsrc/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/',
                     '{hap_reads}.{sequence}.hifiasm' + '.t{}.rsrc'.format(config['num_cpu_high']))
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((12288 + 12288 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 12288 + 12288 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        prefix = lambda wildcards, output: output.hap1_contigs.rsplit('.', 4)[0],
    shell:
        'hifiasm -o {params.prefix} -t {threads} -3 {input.tags_h1} -4 {input.tags_h2} '
            '{input.fastq_h1} {input.fastq_h2} {input.fastq_un} &> {log.hifiasm}'


rule convert_cluster_gfa_to_fasta:
    input:
        hap1 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}/{hap_reads}.{sequence}.dip.hap1.p_ctg.gfa',
        hap2 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}/{hap_reads}.{sequence}.dip.hap2.p_ctg.gfa',
    output:
        fasta_h1 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-hifiasm.h1-un.{sequence}.fasta',
        map_h1 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-hifiasm.h1-un.{sequence}.read-contig.map',
        stats_h1 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-hifiasm.h1-un.{sequence}.contig.stats',
        fasta_h2 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-hifiasm.h2-un.{sequence}.fasta',
        map_h2 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-hifiasm.h2-un.{sequence}.read-contig.map',
        stats_h2 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-hifiasm.h2-un.{sequence}.contig.stats',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-hifiasm.{sequence}.gfa-convert.log'
    benchmark:
        'rsrc/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-hifiasm.{sequence}' + '.t{}.rsrc'.format(config['num_cpu_low'])
    conda:
        '../environment/conda/conda_pyscript.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = lambda wildcards: find_script_path('gfa_to_fasta.py')
    shell:
        '{params.script_exec} --gfa {input.hap1} --n-cpus {threads} '
        '--out-fasta {output.fasta_h1} --out-map {output.map_h1} --out-stats {output.stats_h1} &> {log}'
        ' && '
        '{params.script_exec} --gfa {input.hap2} --n-cpus {threads} '
        '--out-fasta {output.fasta_h2} --out-map {output.map_h2} --out-stats {output.stats_h2} &>> {log}'
