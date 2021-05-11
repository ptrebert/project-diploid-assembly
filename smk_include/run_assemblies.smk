
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


def load_shasta_minreadlength(stats_file, target_cov):
    """
    :param stats_file:
    :param target_cov:
    :return:
    """
    import sys
    min_read_length = 0
    last_delta = sys.maxsize
    with open(stats_file, 'r') as table:
        for line in table:
            if not line.startswith('cov_geq'):
                continue
            rlen_info, cov = line.strip().split()
            rlen = int(rlen_info.split('_')[-1])
            cov = float(cov)
            # note that first entry is geq_0, i.e., no filtering
            # in case the target coverage is not attainable at all
            delta = abs(round(cov - target_cov, 1))
            if delta < last_delta:
                last_delta = delta
                min_read_length = rlen
    return min_read_length


rule derive_shasta_parameter_preset:
    input:
        '{folder_path}/{file_name}.fastq',
        '{folder_path}/{file_name}.stats'
    output:
        '{folder_path}/{file_name}.preset.shasta'
    run:
        import os
        preset = None
        file_name = os.path.basename(input[0])
        _, _, platform_spec = file_name.split('_')[:3]
        target_coverage = int(config['shasta_target_coverage'])
        if platform_spec.startswith('pb'):
            if 'ccs' in platform_spec:
                # https://github.com/chanzuckerberg/shasta/blob/master/conf/PacBio-CCS-Dec2019.conf
                preset = [
                    '[Kmers]',
                    'k = 15',
                    'probability = 0.02',
                    '[MinHash]',
                    'm = 12',
                    'minBucketSize = 20',
                    'maxBucketSize = 100',
                    'minHashIterationCount = 25',
                    'minFrequency = 10',
                    '[ReadGraph]',
                    'maxAlignmentCount = 20',
                    '[Assembly]',
                    'consensusCaller = Modal'
                ]
            elif 'clr' in platform_spec:
                # https://github.com/chanzuckerberg/shasta/blob/master/conf/PacBio-CLR-Dec2019.conf
                preset = [
                    '[ReadGraph]',
                    'maxAlignmentCount = 20',
                    '[Assembly]',
                    'consensusCaller = Modal'
                ]
            else:
                raise ValueError('Cannot handle PacBio platform spec: {} / {}'.format(file_name, platform_spec))
        elif platform_spec.startswith('ont'):
            if platform_spec.endswith('-ul'):
                # https://github.com/chanzuckerberg/shasta/blob/master/conf/Nanopore-UL-Dec2019.conf
                preset = [
                    '[MinHash]',
                    'minBucketSize = 5',
                    'maxBucketSize = 40',
                    'minFrequency = 10',
                    '[Align]',
                    'maxSkip = 60',
                    'maxDrift = 60',
                    'minAlignedMarkerCount = 400',
                    '[Assembly]',
                    'consensusCaller = Bayesian:guppy-3.0.5-a'
                ]
            else:
                # https://github.com/chanzuckerberg/shasta/blob/master/conf/Nanopore-Dec2019.conf
                preset = [
                    '[MinHash]',
                    'minBucketSize = 5',
                    'maxBucketSize = 30',
                    'minFrequency = 5',
                    '[Align]',
                    'minAlignedFraction = 0.4',
                    '[Assembly]',
                    'consensusCaller = Bayesian:guppy-3.0.5-a',
                ]
        else:
            raise ValueError('Unexpected platform spec: {} / {}'.format(file_name, platform_spec))
        min_read_length = load_shasta_minreadlength(input[1], target_coverage)
        preset.append('[Reads]')
        preset.append('minReadLength = ' + str(min_read_length))

        with open(output[0], 'w') as shasta_config:
            _ = shasta_config.write('\n'.join(preset))
    # end of run block


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
        'run/output/reference_assembly/non-hap-res/layout/wtdbg2/{{sample}}_nhr-wtdbg.layout.t{}.rsrc'.format(config['num_cpu_max'])
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
        'run/output/reference_assembly/non-hap-res/{{sample}}_nhr-wtdbg.consensus.t{}.rsrc'.format(config['num_cpu_high'])
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
        'run/output/reference_assembly/non-hap-res/{{sample}}_nhr-flye.layout.t{}.rsrc'.format(config['num_cpu_max'])
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


rule write_peregrine_input_fofn:
    input:
        '{folder_path}/{file_name}.fastq.gz'
    output:
        '{folder_path}/{file_name}.input.peregrine'
    params:
        mount_point = '/wd'
    run:
        with open(output[0], 'w') as fofn:
            _ = fofn.write(params.mount_point + '/' + input[0])


rule compute_peregrine_nonhapres_assembly:
    """
    Why this construct: (yes yes || true)
    1) Peregrine is PacBio software, needs an explicit "yes" to the question
        about the usage scenario (comm. vs. academic/non-comm.)
    2) the command "yes" continues sending a "yes" into the singularity container for all eternity
    3) as soon as Peregrine is done and exits, this leads to a broken pipe
    4) a broken pipe (SIGPIPE) leads to a non-zero exit code for the whole command (specifically, 141)
    5) since Snakemake does not report the exact exit code, one may waste a couple of hours looking
        in the wrong direction for the source of the error...
    6) the " || true" ensures that "yes" also returns 0 as exit code
    7) this is safe-ish since errors from Peregrine (or cp) still signal a failure as desired
    """
    input:
        container = 'output/container/docker/cschin/peregrine_{}.sif'.format(config['peregrine_version']),
        fofn = 'input/fastq/{sample}.input.peregrine'
    output:
        dir_seqdb = directory('output/reference_assembly/non-hap-res/layout/peregrine/{sample}/0-seqdb'),
        dir_index = directory('output/reference_assembly/non-hap-res/layout/peregrine/{sample}/1-index'),
        dir_ovlp = directory('output/reference_assembly/non-hap-res/layout/peregrine/{sample}/2-ovlp'),
        dir_asm = directory('output/reference_assembly/non-hap-res/layout/peregrine/{sample}/3-asm'),
        dir_cns = directory('output/reference_assembly/non-hap-res/layout/peregrine/{sample}/4-cns'),
        assm = 'output/reference_assembly/non-hap-res/{sample}_nhr-pereg.fasta'
    log:
        pereg = 'log/output/reference_assembly/non-hap-res/{sample}_nhr-pereg.log',
        assm_copy = 'log/output/reference_assembly/non-hap-res/{sample}_nhr-pereg.copy.log'
    benchmark:
        os.path.join('run/output/reference_assembly/non-hap-res',
                     '{sample}_nhr-pereg' + '.t{}.rsrc'.format(config['num_cpu_max']))
    envmodules:
        config['env_module_singularity']
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((1499136 if attempt <= 1 else 2949120) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards, attempt: 1499136 if attempt <= 1 else 2949120,
        runtime_hrs = lambda wildcards, attempt: 12 if attempt <= 1 else 16 * attempt
    params:
        bind_folder = lambda wildcards: os.getcwd(),
        out_folder = lambda wildcards, output: os.path.split(output.dir_seqdb)[0],
        singularity = '' if not config.get('env_module_singularity', False) else 'module load {} ; '.format(config['env_module_singularity'])
    shell:
        '{params.singularity}'
        '(yes yes || true) | singularity run --bind {params.bind_folder}:/wd {input.container} '
            ' asm {input.fofn} {threads} {threads} {threads} {threads} {threads} {threads} {threads} {threads} {threads} '
            ' --with-consensus --shimmer-r 3 --best_n_ovlp 8 --output /wd/{params.out_folder} &> {log.pereg} '
            ' && '
            ' cp {output.dir_cns}/cns-merge/ctg_cns.fa {output.assm} &> {log.assm_copy}'


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
        os.path.join('run/output/reference_assembly/non-hap-res',
                     '{sample}_nhr-hifiasm' + '.t{}.rsrc'.format(config['num_cpu_high']))
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((188416 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 188416 * attempt,
        runtime_hrs = lambda wildcards, attempt: 36 * attempt
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
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
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
        fasta = 'output/reference_assembly/non-hap-res/{sample}_nhr-hifiasm.fasta',
        rc_map = 'output/reference_assembly/non-hap-res/{sample}_nhr-hifiasm.read-contig.map',
        stats = 'output/reference_assembly/non-hap-res/{sample}_nhr-hifiasm.contig.stats',
    log:
        'log/output/reference_assembly/non-hap-res/{sample}_nhr-hifiasm.gfa-convert.log'
    benchmark:
        'run/output/reference_assembly/non-hap-res/{sample}_nhr-hifiasm.gfa-convert' + '.t{}.rsrc'.format(config['num_cpu_low'])
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


rule compute_shasta_nonhapres_assembly:
    """
    Non-haplotype resolved assembly = replaces known reference (such as hg38)
    Note that for Shasta, the output directory MUST NOT exist at invocation time, which
    is incompatible with Snakemake's default behavior - rm output dir before start...
    """
    input:
        shasta_exec = 'output/check_files/environment/shasta_version.ok',
        fasta = 'input/fastq/{sample}.fastq',
        config = 'input/fastq/{sample}.preset.shasta'
    output:
        assm_files = expand('output/reference_assembly/non-hap-res/layout/shasta/{{sample}}/Assembly{assm_files}',
                            assm_files=['-BothStrands.gfa', '.gfa', 'GraphChainLengthHistogram.csv',
                                        'Graph-Final.dot', 'Summary.csv', 'Summary.html', 'Summary.json']
                            ),
        aux_files = expand('output/reference_assembly/non-hap-res/layout/shasta/{{sample}}/{aux_files}',
                            aux_files=['Binned-ReadLengthHistogram.csv', 'MarkerGraphEdgeCoverageHistogram.csv',
                                        'MarkerGraphVertexCoverageHistogram.csv', 'PalindromicReads.csv',
                                        'ReadGraphComponents.csv', 'ReadLengthHistogram.csv', 'ReadSummary.csv',
                                        'shasta.conf']
                           ),
        assm_source = 'output/reference_assembly/non-hap-res/layout/shasta/{sample}/Assembly.fasta',
        assembly = 'output/reference_assembly/non-hap-res/{sample}_nhr-shasta.fasta'
    log:
        'log/output/reference_assembly/non-hap-res/{sample}_nhr-shasta.layout.log',
    benchmark:
        'run/output/reference_assembly/non-hap-res/{{sample}}_nhr-shasta.layout.t{}.rsrc'.format(config['num_cpu_max'])
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((1499136 if attempt <= 1 else 2949120) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards, attempt: 1499136 if attempt <= 1 else 2949120,
        runtime_hrs = lambda wildcards, attempt: 12 if attempt <= 1 else 71
    params:
        out_prefix = lambda wildcards, output: os.path.dirname(output.assm_source)
    shell:
        'rm -rfd {params.out_prefix} && '
        'shasta --input {input.fasta} --assemblyDirectory {params.out_prefix} --command assemble '
            ' --memoryMode anonymous --memoryBacking 4K --threads {threads} '
            ' --config {input.config} &> {log}'
            ' && '
            'cp {output.assm_source} {output.assembly}'


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
        'run/output/diploid_assembly/{{variant}}/{{folder_path}}/draft/haploid_assembly/{{hap_reads}}.{{hap}}.wtdbg-consensus.t{}.rsrc'.format(config['num_cpu_high'])
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
        'run/output/diploid_assembly/{{variant}}/{{folder_path}}/draft/haploid_assembly/{{hap_reads}}.{{hap}}.flye.t{}.rsrc'.format(config['num_cpu_max'])
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


rule compute_shasta_haploid_assembly:
    """
    Note that for Shasta, the output directory MUST NOT exist at invocation time, which
    is incompatible with Snakemake's default behavior - rm output dir before start...
    """
    input:
        shasta_exec = 'output/check_files/environment/shasta_version.ok',
        fasta = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_fastq/{hap_reads}.{hap}.fastq',
        config = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_fastq/{hap_reads}.{hap}.preset.shasta',
    output:
        assm_files = expand('output/' + PATH_STRANDSEQ_DGA_JOINT_PROTECTED + '/draft/temp/layout/shasta/{{hap_reads}}.{{hap}}/Assembly{assm_files}',
                            assm_files=['-BothStrands.gfa', '.gfa', 'GraphChainLengthHistogram.csv',
                                        'Graph-Final.dot', 'Summary.csv', 'Summary.html', 'Summary.json']
                            ),
        aux_files = expand('output/' + PATH_STRANDSEQ_DGA_JOINT_PROTECTED + '/draft/temp/layout/shasta/{{hap_reads}}.{{hap}}/{aux_files}',
                            aux_files=['Binned-ReadLengthHistogram.csv', 'MarkerGraphEdgeCoverageHistogram.csv',
                                        'MarkerGraphVertexCoverageHistogram.csv', 'PalindromicReads.csv',
                                        'ReadGraphComponents.csv', 'ReadLengthHistogram.csv', 'ReadSummary.csv',
                                        'shasta.conf']
                           ),
        assm_source = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/temp/layout/shasta/{hap_reads}.{hap}/Assembly.fasta',
        assembly = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_assembly/{hap_reads}-shasta.{hap}.fasta',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_assembly/{hap_reads}-shasta.{hap}.log',
    benchmark:
        os.path.join('run/output', PATH_STRANDSEQ_DGA_JOINT, 'draft/haploid_assembly',
                     '{hap_reads}-shasta.{hap}' + '.t{}.rsrc'.format(config['num_cpu_max']))
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_max']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((1499136 if attempt <= 1 else 2949120) / config['num_cpu_max']),
        mem_total_mb = lambda wildcards, attempt: 1499136 if attempt <= 1 else 2949120,
        runtime_hrs = lambda wildcards, attempt: 12 if attempt <= 1 else 71
    params:
        out_prefix = lambda wildcards, output: os.path.dirname(output.assm_source)
    shell:
        'rm -rfd {params.out_prefix} && '
        'shasta --input {input.fasta} --assemblyDirectory {params.out_prefix} --command assemble '
            ' --memoryMode anonymous --memoryBacking 4K --threads {threads} '
            ' --config {input.config} &> {log}'
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
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/wtdbg2/{{hap_reads}}.{{hap}}.{{sequence}}.wtdbg-layout.t{}.rsrc'.format(config['num_cpu_high'])
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
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{{hap_reads}}.{{hap}}.{{sequence}}.wtdbg-consensus.t{}.rsrc'.format(config['num_cpu_high'])
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
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{{hap_reads}}.{{hap}}.{{sequence}}.flye.t{}.rsrc'.format(config['num_cpu_high'])
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


rule compute_peregrine_haploid_split_assembly:
    """
    Why this construct: (yes yes || true)
    1) Peregrine is PacBio software, needs an explicit "yes" to the question
        about the usage scenario (comm. vs. academic/non-comm.)
    2) the command "yes" continues sending a "yes" into the singularity container for all eternity
    3) as soon as Peregrine is done and exits, this leads to a broken pipe
    4) a broken pipe (SIGPIPE) leads to a non-zero exit code for the whole command (specifically, 141)
    5) since Snakemake does not report the exact exit code, one may waste a couple of hours looking
        in the wrong direction for the source of the error...
    6) the " || true" ensures that "yes" also returns 0 as exit code
    7) this is safe-ish since errors from Peregrine (or cp) still signal a failure as desired
    """
    input:
        container = 'output/container/docker/cschin/peregrine_{}.sif'.format(config['peregrine_version']),
        fofn = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.input.peregrine',
    output:
        dir_seqdb = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/peregrine/{hap_reads}.{hap}.{sequence}/0-seqdb'),
        dir_index = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/peregrine/{hap_reads}.{hap}.{sequence}/1-index'),
        dir_ovlp = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/peregrine/{hap_reads}.{hap}.{sequence}/2-ovlp'),
        dir_asm = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/peregrine/{hap_reads}.{hap}.{sequence}/3-asm'),
        dir_cns = directory('output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/peregrine/{hap_reads}.{hap}.{sequence}/4-cns'),
        assm = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-pereg.{hap}.{sequence}.fasta',
    log:
        pereg = 'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-pereg.{hap}.{sequence}.log',
        assm_copy = 'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-pereg.{hap}.{sequence}.copy.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{{hap_reads}}-pereg.{{hap}}.{{sequence}}.t{}.rsrc'.format(config['num_cpu_high'])
    envmodules:
        config['env_module_singularity']
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((110592 if attempt <= 1 else 188416 * (attempt - 1)) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 110592 if attempt <= 1 else 188416 * (attempt - 1),
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        bind_folder = lambda wildcards: os.getcwd(),
        out_folder = lambda wildcards, output: os.path.split(output.dir_seqdb)[0],
        singularity = '' if not config.get('env_module_singularity', False) else 'module load {} ; '.format(config['env_module_singularity'])
    shell:
        '{params.singularity}'
        '(yes yes || true) | singularity run --bind {params.bind_folder}:/wd {input.container} '
            ' asm {input.fofn} {threads} {threads} {threads} {threads} {threads} {threads} {threads} {threads} {threads} '
            ' --with-consensus --shimmer-r 3 --best_n_ovlp 8 --output /wd/{params.out_folder} &> {log.pereg} '
            ' && '
            ' cp {output.dir_cns}/cns-merge/ctg_cns.fa {output.assm} &> {log.assm_copy}'


rule compute_shasta_haploid_split_assembly:
    """
    Note that for Shasta, the output directory MUST NOT exist at invocation time, which
    is incompatible with Snakemake's default behavior - rm output dir before start...
    """
    input:
        shasta_exec = 'output/check_files/environment/shasta_version.ok',
        fasta = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq',
        config = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.preset.shasta',
    output:
        assm_files = expand('output/' + PATH_STRANDSEQ_DGA_SPLIT_PROTECTED + '/draft/temp/layout/shasta/{{hap_reads}}.{{hap}}.{{sequence}}/Assembly{assm_files}',
                            assm_files=['-BothStrands.gfa', '.gfa', 'GraphChainLengthHistogram.csv',
                                        'Graph-Final.dot', 'Summary.csv', 'Summary.html', 'Summary.json']
                            ),
        aux_files = expand('output/' + PATH_STRANDSEQ_DGA_SPLIT_PROTECTED + '/draft/temp/layout/shasta/{{hap_reads}}.{{hap}}.{{sequence}}/{aux_files}',
                            aux_files=['Binned-ReadLengthHistogram.csv', 'MarkerGraphEdgeCoverageHistogram.csv',
                                        'MarkerGraphVertexCoverageHistogram.csv', 'PalindromicReads.csv',
                                        'ReadGraphComponents.csv', 'ReadLengthHistogram.csv', 'ReadSummary.csv',
                                        'shasta.conf']
                           ),
        assm_source = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/shasta/{hap_reads}.{hap}.{sequence}/Assembly.fasta',
        assembly = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-shasta.{hap}.{sequence}.fasta',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-shasta.{hap}.{sequence}.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{{hap_reads}}-shasta.{{hap}}.{{sequence}}.t{}.rsrc'.format(config['num_cpu_high'])
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((49152 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        out_prefix = lambda wildcards, output: os.path.dirname(output.assm_source)
    shell:
        'rm -rfd {params.out_prefix} && '
        'shasta --input {input.fasta} --assemblyDirectory {params.out_prefix} --command assemble '
            ' --memoryMode anonymous --memoryBacking 4K --threads {threads} '
            ' --config {input.config} &> {log}'
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
        os.path.join('run/output', PATH_STRANDSEQ_DGA_SPLIT, 'draft/haploid_assembly/',
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
        hap1_contigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}/{hap_reads}.{sequence}.hap1.p_ctg.gfa',
        hap2_contigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}/{hap_reads}.{sequence}.hap2.p_ctg.gfa',
        raw_unitigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}/{hap_reads}.{sequence}.dip.r_utg.gfa',
        discard = multiext(
                'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}/{hap_reads}.{sequence}',
                '.hap1.p_ctg.noseq.gfa', '.hap2.p_ctg.noseq.gfa',
                '.dip.r_utg.noseq.gfa'
            )
        # as above: exclude read self overlaps from output
        # '.ec.bin', '.ovlp.reverse.bin', '.ovlp.source.bin',
    log:
        hifiasm = 'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}.hifiasm.log',
    benchmark:
        os.path.join('run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/',
                     '{hap_reads}.{sequence}.hifiasm' + '.t{}.rsrc'.format(config['num_cpu_high']))
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((12288 + 12288 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 12288 + 12288 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        prefix = lambda wildcards, output: output.hap1_contigs.rsplit('.', 3)[0],
    shell:
        'hifiasm -o {params.prefix} -t {threads} -3 {input.tags_h1} -4 {input.tags_h2} '
            '{input.fastq_h1} {input.fastq_h2} {input.fastq_un} &> {log.hifiasm}'


rule convert_cluster_gfa_to_fasta:
    input:
        hap1 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}/{hap_reads}.{sequence}.hap1.p_ctg.gfa',
        hap2 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/temp/layout/hifiasm/{hap_reads}.{sequence}/{hap_reads}.{sequence}.hap2.p_ctg.gfa',
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
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-hifiasm.{sequence}' + '.t{}.rsrc'.format(config['num_cpu_low'])
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
