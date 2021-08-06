
localrules: dump_reads_fofn

# raw Illumina
# ec Illumina
# raw HiFi
# af HiFi
# ec HiFi
# raw ONT
# ONTec w/ single parameterization?

INPUT_READS = [
        'input/ont/NA18989_ONTUL_guppy-5.0.11-sup-prom.fasta.gz',
        'input/ont/NA18989_ONTUL_guppy-4.0.11-hac-prom.fasta.gz',
        'input/hifi/NA18989_HIFIEC_hifiasm-v0.15.4.fasta.gz',
        'input/hifi/NA18989_HIFIAF_pgas-v14-dev.fastq.gz',
        'input/short/NA18989_ERR3239679.fasta.gz'
    ]


READSETS = [os.path.basename(x.rsplit('.', 2)[0]) for x in INPUT_READS]

wildcard_constraints:
    sample = 'NA18989'


def find_script_path(script_name, subfolder=''):
    import os

    current_root = workflow.basedir
    last_root = ''

    script_path = None

    for _ in range(workflow.basedir.count('/')):
        if last_root.endswith('project-diploid-assembly'):
            raise RuntimeError('Leaving project directory tree (next: {}). '
                               'Cannot find script {} (subfolder: {}).'.format(current_root, script_name, subfolder))
        check_path = os.path.join(current_root, 'scripts', subfolder).rstrip('/')  # if subfolder is empty string
        if os.path.isdir(check_path):
            check_script = os.path.join(check_path, script_name)
            if os.path.isfile(check_script):
                script_path = check_script
                break
        last_root = current_root
        current_root = os.path.split(current_root)[0]

    if script_path is None:
        raise RuntimeError('Could not find script {} (subfolder {}). '
                           'Started at path: {}'.format(script_name, subfolder, workflow.basedir))
    return script_path


rule run_all:
    input:
        readsets = INPUT_READS,
        ontul_stats = expand(
            'output/read_stats/input/{sample}_ONTUL_{basecaller}.summary.tsv',
            sample=['NA18989'],
            basecaller=['guppy-5.0.11-sup-prom', 'guppy-4.0.11-hac-prom']
        ),
        kmer_dbs = expand('output/kmer_db/{readset}.total.count', readset=READSETS),
        ont_align = expand(
            'output/alignments/ont_to_mbg_hifi/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.gaf',
            zip,
            sample=['NA18989'] * 4,
            readset=['ONTUL_guppy-5.0.11-sup-prom', 'ONTUL_guppy-5.0.11-sup-prom', 'ONTUL_guppy-4.0.11-hac-prom', 'ONTUL_guppy-4.0.11-hac-prom'],
            kmer=[1001, 3001] * 2,
            window=[500, 2000] * 2
        ),
        ukm5sup = 'output/kmer_stats/NA18989_ONTUL_guppy-5.0.11-sup-prom.unsupported.counts.tsv',
        ukm4hac = 'output/kmer_stats/NA18989_ONTUL_guppy-4.0.11-hac-prom.unsupported.counts.tsv',
        align_stats = expand(
            'output/aln_stats/ont_to_mbg_hifi/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.gaf-stats.tsv.gz',
            zip,
            sample=['NA18989'] * 4,
            readset=['ONTUL_guppy-5.0.11-sup-prom', 'ONTUL_guppy-5.0.11-sup-prom', 'ONTUL_guppy-4.0.11-hac-prom', 'ONTUL_guppy-4.0.11-hac-prom'],
            kmer=[1001, 3001] * 2,
            window=[500, 2000] * 2
        ),
        ontec_stats = expand(
            'output/read_stats/ontec/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.clip-ec.summary.tsv',
            zip,
            sample=['NA18989'] * 4,
            readset=['ONTUL_guppy-5.0.11-sup-prom', 'ONTUL_guppy-5.0.11-sup-prom', 'ONTUL_guppy-4.0.11-hac-prom', 'ONTUL_guppy-4.0.11-hac-prom'],
            kmer=[1001, 3001] * 2,
            window=[500, 2000] * 2
        ),
        ontec_ukm = expand(
            'output/kmer_stats/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.unsupported.counts.tsv',
            zip,
            sample=['NA18989'] * 4,
            readset=['ONTUL_guppy-5.0.11-sup-prom', 'ONTUL_guppy-5.0.11-sup-prom', 'ONTUL_guppy-4.0.11-hac-prom', 'ONTUL_guppy-4.0.11-hac-prom'],
            kmer=[1001, 3001] * 2,
            window=[500, 2000] * 2
        ),
        ontul_jaccard = expand(
            'output/kmer_stats/{sample}_ONTUL_{basecaller1}_vs_{sample}_ONTUL_{basecaller2}.union.count',
            sample=['NA18989'],
            basecaller1=['guppy-5.0.11-sup-prom'],
            basecaller2=['guppy-4.0.11-hac-prom']
        ),
        ontec_jaccard = expand(
            'output/kmer_stats/{readset1}_vs_{readset2}.union.count',
            zip,
            readset1=[
                'NA18989_ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k1001-w500',
                'NA18989_HIFIEC_hifiasm-v0.15.4',
                'NA18989_ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k1001-w500',
                'NA18989_ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k3001-w2000',
                'NA18989_ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k1001-w500',
                'NA18989_ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k3001-w2000',
                'NA18989_HIFIEC_hifiasm-v0.15.4',
            ],
            readset2=[
                'NA18989_ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k3001-w2000',
                'NA18989_HIFIAF_pgas-v14-dev',
                'NA18989_HIFIEC_hifiasm-v0.15.4',
                'NA18989_HIFIEC_hifiasm-v0.15.4',
                'NA18989_ONTUL_guppy-4.0.11-hac-prom_MAP-TO_HIFIEC.mbg-k1001-w500',
                'NA18989_ONTUL_guppy-4.0.11-hac-prom_MAP-TO_HIFIEC.mbg-k3001-w2000',
                'NA18989_ERR3239679'
            ]
        )


def select_ont_input(wildcards):

    if wildcards.basecaller == 'guppy-5.0.11-sup-prom':
        fastq = [
            '/gpfs/project/projects/medbioinf/data/hgsvc_ontul/NA18989/20210510_210428_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom_fastq_pass.fastq.gz',
            '/gpfs/project/projects/medbioinf/data/hgsvc_ontul/NA18989/20210519_210512_21-lee-006_PCT0053_2-A9-D9_guppy-5.0.11-sup-prom_fastq_pass.fastq.gz'
        ]
    elif wildcards.basecaller == 'guppy-4.0.11-hac-prom':
        fastq = [
            '/gpfs/project/projects/medbioinf/data/hgsvc_ontul/NA18989_JAX_old-basecaller/GM18989_GT21-08680.20210510_guppy-4.0.11-hac-prom_fastq_pass.fastq.gz',
            '/gpfs/project/projects/medbioinf/data/hgsvc_ontul/NA18989_JAX_old-basecaller/GM18989_GT21-08680.20210519_guppy-4.0.11-hac-prom_fastq_pass.fastq.gz'
        ]
    else:
        raise ValueError(str(wildcards))
    return fastq


rule merge_ont_data:
    input:
        select_ont_input
    output:
        'input/ont/NA18989_ONTUL_{basecaller}.fasta.gz'
    benchmark:
        'rsrc/input/ont/NA18989_ONTUL_{basecaller}.merge.rsrc'
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    shell:
        'pigz -d -c {input} | seqtk seq -A -C | pigz -p {threads} --best > {output}'


rule compute_input_read_stats:
    input:
        'input/ont/{sample}_ONTUL_{basecaller}.fasta.gz'
    output:
        dump = 'output/read_stats/input/{sample}_ONTUL_{basecaller}.dump.pck',
        summary = 'output/read_stats/input/{sample}_ONTUL_{basecaller}.summary.tsv',
    log:
        'log/output/read_stats/{sample}_ONTUL_{basecaller}.comp-stats.log'
    benchmark:
        'rsrc/output/read_stats/{sample}_ONTUL_{basecaller}.comp-stats.rsrc'
    conda: '../../../environment/conda/conda_pyscript.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt * attempt
    params:
        script_exec = find_script_path('collect_read_stats.py')
    shell:
        '{params.script_exec} --input-files {input} --output {output.dump} --summary-output {output.summary} '
            '--genome-size 3100000000 --num-cpu {threads} &> {log}'


rule dump_reads_fofn:
    input:
        INPUT_READS
    output:
        'input/NA18989_reads.fofn'
    run:
        with open(output[0], 'w') as dump:
            for record in INPUT_READS:
                _ = dump.write(record + '\n')


rule build_bifrost_colored_dbg:
    input:
        container = 'bifrost_x86-64_AVX2_k64.sif',
        read_fofn = 'input/NA18989_reads.fofn'
    output:
        'output/cdbg/NA18989.gfa',
        'output/cdbg/NA18989.bfg_colors'
    log:
       'log/output/cdbg/NA18989.build.log',
    benchmark:
        'rsrc/output/cdbg/NA18989.build.rsrc',
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 180224 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    params:
        kmer_size = 63,
        out_prefix = lambda wildcards, output: output[0].rsplit('.', 1)[0]
    shell:
        'module load Singularity && singularity exec {input.container} '
        'Bifrost build --input-seq-file {input.read_fofn} '
        '--output-file {params.out_prefix} --threads {threads} --colors --kmer-length {params.kmer_size} '
        '--verbose &> {log}'


def count_kmer_runtime(wildcards, attempt):

    if 'HIFI' in wildcards.readset:
        return 24 * attempt
    elif 'ONT' in wildcards.readset:
        return 12 * attempt
    else:
        return attempt * attempt * attempt


def count_kmer_memory(wildcards, attempt, unit='mb'):

    if 'HIFI' in wildcards.readset:
        mem = 176128
    elif 'ONT' in wildcards.readset:
        mem = 90112
    else:
        mem = 32768
    if unit == 'gb':
        mem = int(mem / 1024)
    return mem * attempt


def select_sequence_input(wildcards):

    readset = [x for x in INPUT_READS if wildcards.readset in x]
    if len(readset) == 0:
        if 'mbg' in wildcards.readset and 'HIFIEC' in wildcards.readset:
            ontec_path = 'output/alignments/ont_to_mbg_hifi/' + wildcards.readset + '.clip-ec.fa.gz'
            readset = [ontec_path]
        else:
            raise ValueError(f'Cannot determine desired output: {str(wildcards)}')
    elif len(readset) == 1:
        pass
    else:
        raise ValueError(f'Cannot determine desired output: {str(wildcards)}')
    return readset[0]


rule meryl_count_kmers:
    input:
        sequence = select_sequence_input
    output:
        kmer_db = directory('output/kmer_db/{readset}.meryl'),
        total = 'output/kmer_db/{readset}.total.count',
        distinct = 'output/kmer_db/{readset}.distinct.count',
        singleton = 'output/kmer_db/{readset}.singleton.count',
    benchmark:
        'rsrc/output/kmer_db/{readset}.meryl.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: count_kmer_memory(wildcards, attempt),
        mem_total_gb = lambda wildcards, attempt: count_kmer_memory(wildcards, attempt, 'gb'),
        runtime_hrs = lambda wildcards, attempt: count_kmer_runtime(wildcards, attempt)
    params:
        kmer_size = 31,
        use_hpc = 'compress'
    shell:
        "meryl count k={params.kmer_size} threads={threads} memory={resources.mem_total_gb} {params.use_hpc} output {output.kmer_db} {input.sequence} && "
            "meryl print {output.kmer_db} | cut -f 2 | awk '{{s+=$1}} END {{print s}}' > {output.total} && "
            "meryl print {output.kmer_db} | wc -l > {output.distinct} && "
            "meryl print [equal-to 1 {output.kmer_db}] | wc -l > {output.singleton} "


rule build_hifi_read_dbg:
    """
    MBG v1.05+ has better mem management, should limit <100G
    """
    input:
        'input/hifi/NA18989_HIFIEC_hifiasm-v0.15.4.fasta.gz',
    output:
        'output/mbg_hifi/NA18989_HIFIEC.mbg-k{kmer}-w{window}.gfa'
    log:
        'log/output/mbg_hifi/NA18989_HIFIEC.mbg-k{kmer}-w{window}.log'
    benchmark:
        'rsrc/output/mbg_hifi/NA18989_HIFIEC.mbg-k{kmer}-w{window}.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 65536 + 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 8 * attempt
    shell:
        'MBG -i {input} -o {output} -t {threads} -k {wildcards.kmer} -w {wildcards.window} &> {log}'


rule ont_to_graph_alignment:
    """
    """
    input:
        container = ancient('graphaligner.sif'),
        graph = 'output/mbg_hifi/{sample}_HIFIEC.mbg-k{kmer}-w{window}.gfa',
        reads = 'input/ont/{sample}_{readset}.fasta.gz',
    output:
        gaf = 'output/alignments/ont_to_mbg_hifi/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.gaf',
        ec_reads_clip = 'output/alignments/ont_to_mbg_hifi/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.clip-ec.fa.gz',
    log:
        'log/output/alignments/ont_to_mbg_hifi/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.ga.log'
    benchmark:
        'rsrc/output/alignments/ont_to_mbg_hifi/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.ga.rsrc'
#    conda:
#        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 90112 + 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: 36 * attempt
    params:
        preset = 'dbg'
    shell:
        'module load Singularity && singularity exec {input.container} '
        'GraphAligner -g {input.graph} -f {input.reads} '
            '-x {params.preset} -t {threads} '
            '--min-alignment-score 10000 --multimap-score-fraction 1 '
            '--corrected-clipped-out {output.ec_reads_clip} '
            '-a {output.gaf} &> {log}'


rule compute_ontec_read_stats:
    input:
        'output/alignments/ont_to_mbg_hifi/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.clip-ec.fa.gz',
    output:
        dump = 'output/read_stats/ontec/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.clip-ec.dump.pck',
        summary = 'output/read_stats/ontec/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.clip-ec.summary.tsv',
    log:
        'log/output/read_stats/ontec/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.clip-ec.comp-stats.log'
    benchmark:
        'rsrc/output/read_stats/ontec/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.clip-ec.comp-stats.rsrc'
    conda: '../../../environment/conda/conda_pyscript.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt * attempt
    params:
        script_exec = find_script_path('collect_read_stats.py')
    shell:
        '{params.script_exec} --input-files {input} --output {output.dump} --summary-output {output.summary} '
            '--genome-size 3100000000 --num-cpu {threads} &> {log}'


rule dump_unsupported_kmers:
    """
    Piping the k-mer counts to pigz for direct compression
    extensively prolongs runtime (days...)
    """
    input:
        ont_db = 'output/kmer_db/{readset}.meryl',
        hifiec = 'output/kmer_db/NA18989_HIFIEC_hifiasm-v0.15.4.meryl',
        hifiaf = 'output/kmer_db/NA18989_HIFIAF_pgas-v14-dev.meryl',
        short = 'output/kmer_db/NA18989_ERR3239679.meryl',
    output:
        'output/kmer_stats/{readset}.unsupported.txt'
    benchmark:
        'rsrc/output/kmer_stats/{readset}.unsupported.meryl.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: 16 * attempt
    shell:
        'meryl print [difference {input.ont_db} {input.hifiec} {input.hifiaf} {input.short}] > {output}'


rule create_unsupported_kmer_histogram:
    input:
        'output/kmer_stats/{readset}.unsupported.txt'
    output:
        'output/kmer_stats/{readset}.unsupported.counts.tsv'
    resources:
        runtime_hrs = lambda wildcards, attempt: 16 * attempt
    run:
        import collections as col

        hist = col.Counter()

        with open(input[0], 'r') as dump:
            for ln, line in enumerate(dump, start=1):
                try:
                    _, abundance = line.split()
                    hist[int(abundance)] += 1
                except (ValueError, IndexError) as err:
                    if not line.strip():
                        continue
                    raise ValueError(f'LN: {ln} - {str(err)}')
        
        with open(output[0], 'w') as dump:
            _ = dump.write('abundance\tcount\n')
            for k in sorted(hist.keys()):
                c = hist[k]
                _ = dump.write(f'{k}\t{c}\n')
    # END OF RUN BLOCK


rule prepare_ont_any_jaccard:
    input:
        db1 = 'output/kmer_db/{readset1}.meryl',
        db2 = 'output/kmer_db/{readset2}.meryl',
    output:
        union = 'output/kmer_stats/{readset1}_vs_{readset2}.union.count',
        intersect = 'output/kmer_stats/{readset1}_vs_{readset2}.intersect.count',
    benchmark:
        'rsrc/output/kmer_stats/{readset1}_vs_{readset2}.prepjacc.rsrc',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: 16 * attempt
    shell:
        'meryl print [union {input.db1} {input.db2}] | wc -l > {output.union} '
        ' && '
        'meryl print [intersect {input.db1} {input.db2}] | wc -l > {output.intersect} '


def process_gaf_line(translation_table, gaf_line):

    _, rdlen, alns, alne, _, path, _, _, _, resmatch, _, _, _, ascore, _, identity, _ = gaf_line.split()
    aligned_length = int(alne) - int(alns)
    aligned_fraction = round(aligned_length / int(rdlen) * 100, 2)
    matched_fraction = round(int(resmatch) / int(rdlen) * 100, 2)
    ascore = int(round(float(ascore.split(':')[-1]), 0))
    identity = round(float(identity.split(':')[-1]) * 100, 2)
    nodes = path.translate(translation_table).split()
    return rdlen, aligned_length, aligned_fraction, matched_fraction, ascore, identity, nodes


rule collect_gaf_statistics:
    input:
        'output/alignments/ont_to_mbg_hifi/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.gaf',
    output:
        'output/aln_stats/ont_to_mbg_hifi/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.gaf-stats.tsv.gz',
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt * attempt
    run:
        import gzip

        translation_table = dict((i,i) for i in '1234567890')
        translation_table['>'] = ' '
        translation_table['<'] = ' '
        translation_table = str.maketrans(translation_table)

        data_keys = [
            'read_length',
            'aligned_length',
            'aligned_fraction',
            'matched_fraction',
            'tag_AS',
            'tag_ID',
            'node'
        ]

        with gzip.open(output[0], 'wt') as dump:
            _ = dump.write(f'key\tvalue\n')
            with open(input[0], 'r') as gaf:
                for ln, line in enumerate(gaf, start=1):
                    values = process_gaf_line(translation_table, line)
                    for key, value in zip(data_keys, values):
                        if key == 'node':
                            for node in value:
                                _ = dump.write(f'node\t{node}\n')
                        else:
                            _ = dump.write(f'{key}\t{value}\n')

    # END OF RUN BLOCK
