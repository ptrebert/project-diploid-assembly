
localrules: dump_reads_fofn, run_all, compute_global_query_qv_estimate

ruleorder: meryl_query_only_kmer_db > meryl_count_kmers

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
        kmer_dbs = expand('output/kmer_stats/{readset}.meryl.stats.h5', readset=READSETS),
        ont_align = expand(
            'output/alignments/ont_to_mbg_hifi/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.gaf',
            zip,
            sample=['NA18989'] * 4,
            readset=['ONTUL_guppy-5.0.11-sup-prom', 'ONTUL_guppy-5.0.11-sup-prom', 'ONTUL_guppy-4.0.11-hac-prom', 'ONTUL_guppy-4.0.11-hac-prom'],
            kmer=[1001, 3001] * 2,
            window=[500, 2000] * 2
        ),
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
        ),
        cache_read_cov = expand(
            'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.cache.h5',
            sample=['NA18989'],
            readset=[
                'ONTUL_guppy-5.0.11-sup-prom',
                'HIFIEC_hifiasm-v0.15.4',
                'ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k1001-w500',
                'ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k3001-w2000',
            ],
            reference=['T2Tv11_38p13Y_chm13']
        ),
        global_qv_est = expand(
            'output/qv_estimate/{sample}_{readset}_REF_{short_reads}.qv.tsv',
            sample=['NA18989'],
            readset=[
                'HIFIEC_hifiasm-v0.15.4',
                'HIFIAF_pgas-v14-dev',
                'ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k1001-w500',
                'ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k3001-w2000',
                'ONTUL_guppy-5.0.11-sup-prom'
            ],
            short_reads=['ERR3239679']
        ),
        read_qv_est = expand(
            'output/qv_estimate/{sample}_{readset}_REF_{short_reads}.seq-qv.h5',
            sample=['NA18989'],
            readset=[
                'HIFIEC_hifiasm-v0.15.4',
                'HIFIAF_pgas-v14-dev',
                'ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k1001-w500',
                'ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k3001-w2000',
                'ONTUL_guppy-5.0.11-sup-prom'
            ],
            short_reads=['ERR3239679']
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


rule short_read_quality_trimming:
    """
    Note that due to the rather "hidden" way trim-galore/cutadapt
    are handling multithreading, the number of threads "4" is
    internally translated to something like 15 or 16 (reading, writing etc.)
    """
    input:
        mate1 = 'input/short/{readset}_1.fastq.gz',
        mate2 = 'input/short/{readset}_2.fastq.gz'
    output:
        mate1 = 'input/short/{readset}/trimmed/{readset}_1_val_1.fq.gz',
        report1 = 'input/short/{readset}/trimmed/{readset}_1.fastq.gz_trimming_report.txt',
        mate2 = 'input/short/{readset}/trimmed/{readset}_2_val_2.fq.gz',
        report2 = 'input/short/{readset}/trimmed/{readset}_2.fastq.gz_trimming_report.txt'
    log:
        'log/input/short/{readset}.trimming.log'
    benchmark:
        'run/input/short/{readset}.trimming.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 16
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * 12,
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
    params:
        quality_trim = 20,
        min_read_length = 51,
        outdir = lambda wildcards, input: os.path.join('input', 'short', wildcards.readset, 'trimmed')
    shell:
        'trim_galore --quality {params.quality_trim} --length {params.min_read_length} '
        '--trim-n --output_dir {params.outdir} --cores 4 --paired '
        '{input.mate1} {input.mate2} &> {log}'


def compute_lighter_alpha(trim_report1, trim_report2, genomesize):
    """
    rule of thumb: alpha=(7/C)
    :param trim_report1:
    :param trim_report2:
    :param genomesize:
    :return:
    """
    num_bp1 = 0
    num_bp2 = 0
    if not all([os.path.isfile(x) for x in [trim_report1, trim_report2]]):
        # could be dry run
        return -1
    with open(trim_report1, 'r') as report:
        for line in report:
            if not line.startswith('Total written'):
                continue
            parts = line.strip().split()
            num_bp1 = int(parts[-3].replace(',', ''))
            break

    with open(trim_report2, 'r') as report:
        for line in report:
            if not line.startswith('Total written'):
                continue
            parts = line.strip().split()
            num_bp2 = int(parts[-3].replace(',', ''))
            break
    total_cov = num_bp1 + num_bp2
    avg_cov = total_cov / genomesize
    alpha = round(7/avg_cov, 3)
    return alpha


rule short_read_error_correction:
    input:
        mate1 = 'input/short/{readset}/trimmed/{readset}_1_val_1.fq.gz',
        report1 = 'input/short/{readset}/trimmed/{readset}_1.fastq.gz_trimming_report.txt',
        mate2 = 'input/short/{readset}/trimmed/{readset}_2_val_2.fq.gz',
        report2 = 'input/short/{readset}/trimmed/{readset}_2.fastq.gz_trimming_report.txt'
    output:
        mate1 = 'input/short/{readset}/corrected/{readset}_1_val_1.cor.fq.gz',
        mate2 = 'input/short/{readset}/corrected/{readset}_2_val_2.cor.fq.gz',
    log:
        'log/input/short/{readset}.corr.log'
    benchmark:
        'run/input/short/{readset}.corr.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: 8 * attempt,
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt
    params:
        kmer_size = 31,  # this k-mer size is used for consistency with existing meryl DBs
        alpha = lambda wildcards, input: compute_lighter_alpha(input.report1, input.report2, int(6.2e9)),
        genomesize = int(6.2e9),
        outdir = lambda wildcards, output: os.path.dirname(output.mate1)
    shell:
        'lighter -r {input.mate1} -r {input.mate2} '
        '-k {params.kmer_size} {params.genomesize} {params.alpha} '
        '-od {params.outdir} -t {threads} -zlib 6 &> {log}'


rule merge_error_corrected_short_reads:
    """
    This rule is more convenience s.t. all downstream rules can stay as-is
    """
    input:
        mate1 = 'input/short/{readset}/corrected/{readset}_1_val_1.cor.fq.gz',
        mate2 = 'input/short/{readset}/corrected/{readset}_2_val_2.cor.fq.gz',
    output:
        'input/short/{readset}.fasta.gz'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: 24 * attempt
    shell:
        'pigz -d -c {input} | seqtk seq -A -C | pigz -p {threads} --best > {output}'


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


rule meryl_count_reference_kmers:
    input:
        sequence = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta')
    output:
        kmer_db = directory('/gpfs/project/projects/medbioinf/data/references/{reference}.k15.db'),
        rep_kmer = '/gpfs/project/projects/medbioinf/data/references/{reference}.k15.rep-grt09998.txt'
    benchmark:
        'rsrc/output/kmer_db/{reference}.meryl-ref-kmer.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 12288 + 4096 * attempt,
        mem_total_gb = lambda wildcards, attempt: int((12288 + 4096 * attempt)/1024) - 2,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        "meryl count k=15 threads={threads} memory={resources.mem_total_gb} output {output.kmer_db} {input.sequence} "
            " && "
            "meryl print greater-than distinct=0.9998 {output.kmer_db} > {output.rep_kmer}"


rule meryl_count_kmers:
    input:
        sequence = select_sequence_input
    output:
        kmer_db = directory('output/kmer_db/{readset}.meryl'),
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
        "meryl count k={params.kmer_size} threads={threads} memory={resources.mem_total_gb} {params.use_hpc} output {output.kmer_db} {input.sequence}"


rule meryl_dump_db_statistics:
    """
    This operation is not documented in the cli help
    """
    input:
        kmer_db = 'output/kmer_db/{readset}.meryl'
    output:
        stats = 'output/kmer_stats/{readset}.meryl.statistics'
    benchmark:
        'rsrc/output/kmer_stats/{readset}.meryl.stats.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        "meryl statistics {input.kmer_db} > {output.stats}"


rule store_meryl_db_statistics:
    """
    NB: unique = singletons, k-mers with abundance 1

    FORMAT:

    Found 1 command tree.
    Number of 31-mers that are:
        unique              698268532  (exactly one instance of the kmer is in the input)
        distinct           2560577684  (non-redundant kmer sequences in the input)
        present           70996266663  (...)
        missing   4611686015866810220  (non-redundant kmer sequences not in the input)

                number of   cumulative   cumulative     presence
                distinct     fraction     fraction   in dataset
    frequency        kmers     distinct        total       (1e-6)
    --------- ------------ ------------ ------------ ------------
            1    698268532       0.2727       0.0098     0.000014
            2     20636896       0.2808       0.0104     0.000028
            3      6183975       0.2832       0.0107     0.000042
            4      3214166       0.2844       0.0109     0.000056
            5      2299957       0.2853       0.0110     0.000070
    """
    input:
        stats = 'output/kmer_stats/{readset}.meryl.statistics'
    output:
        hdf = 'output/kmer_stats/{readset}.meryl.stats.h5'
    benchmark:
        'rsrc/output/kmer_stats/{readset}.meryl.hdf.rsrc'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    run:
        import pandas as pd
        db_statistics = ['unique', 'distinct', 'present', 'missing']
        db_stats = {}
        kmer_freqs = []
        freq_line = False

        with open(input.stats, 'r') as txt:
            for ln, line in enumerate(txt, start=1):
                if not line.strip():
                    continue
                elif freq_line:
                    columns = line.strip().split()
                    assert len(columns) == 5, f'Malformed frequency line {ln}: {line.strip()}'
                    kmer_freqs.append(
                        (
                            int(columns[0]),
                            int(columns[1]),
                            float(columns[2]),
                            float(columns[3]),
                            float(columns[4])
                        )
                    )
                elif line.strip().startswith('Number of'):
                    parts = line.strip().split()
                    kmer_size = int(parts[2].split('-')[0])
                    db_stats['kmer_size'] = kmer_size
                else:
                    if any(line.strip().startswith(x) for x in db_statistics):
                        parts = line.strip().split()
                        statistic = parts[0]
                        assert statistic in db_statistics
                        try:
                            value = int(parts[1])
                        except ValueError:
                            # this happens b/c of "distinct" appearing
                            # twice (also in table header)
                            if parts[1] == 'fraction':
                                continue
                            raise
                        db_stats[statistic] = value
                    elif line.strip().startswith('-------'):
                        freq_line = True
                        continue
                    else:  # e.g., table header
                        continue
        assert len(db_stats) == 5, f'Missing DB stats: {db_stats}'
        db_stats = pd.Series(db_stats, index=None, name='db_statistics')
        kmer_freqs = pd.DataFrame.from_records(
            kmer_freqs,
            columns=[
                'frequency',
                'num_distinct_kmers',
                'cum_fraction_distinct',
                'cum_fraction_total',
                'presence_in_dataset'
            ]
        )
        with pd.HDFStore(output.hdf, 'w', complevel=9) as hdf:
            hdf.put('statistics', db_stats, format='fixed')
            hdf.put('kmer_freq', kmer_freqs, format='fixed')
    # END OF RUN BLOCK


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


def select_winnowmap_reads(wildcards):

    if 'MAP-TO_HIFIEC' in wildcards.readset:
        reads = [f'output/alignments/ont_to_mbg_hifi/{wildcards.sample}_{wildcards.readset}.clip-ec.fa.gz']
    elif 'ONTUL' in wildcards.readset:
        reads = [f'input/ont/{wildcards.sample}_{wildcards.readset}.fasta.gz']
    elif 'HIFIEC' in wildcards.readset:
        reads = [f'input/hifi/{wildcards.sample}_{wildcards.readset}.fasta.gz']
    elif 'HIFIAF' in wildcards.readset:
        reads = [f'input/hifi/{wildcards.sample}_{wildcards.readset}.fastq.gz']
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
    """
    aligning the ONT-UL reads is driving me crazy...
    now use the hammer
    """
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


rule ontqc_wmap_align_readsets:
    input:
        reads = select_winnowmap_reads,
        reference = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta'),
        ref_repkmer = '/gpfs/project/projects/medbioinf/data/references/{reference}.k15.rep-grt09998.txt'
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
        preset = lambda wildcards: set_winnowmap_preset(wildcards)
    shell:
        'winnowmap -W {input.ref_repkmer} -k 15 -t {resources.align_threads} -Y -L --eqx --MD -a -x {params.preset} '
        '-R "@RG\\tID:{params.readgroup_id}\\tSM:{params.individual}" --secondary=no '
        '{input.reference} {input.reads} | '
        'samtools view -u -F 260 | '
        'samtools sort -m {resources.mem_sort_mb}M -@ {resources.sort_threads} -O BAM > {output.bam} '
        ' && '
        'samtools index {output.bam}'


rule dump_alignment_coverage:
    input:
        bam = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.psort.bam',
        bai = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.psort.bam.bai',
    output:
        bg = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.cov.bedgraph',
    benchmark:
        'rsrc/output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.cov.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 1
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: 12 * attempt,
    shell:
        'bedtools genomecov -bg -ibam {input.bam} > {output.bg}'


rule cache_read_coverage:
    input:
        bg = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.cov.bedgraph',
        fai = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta.fai'),
    output:
        hdf = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.cache.h5',
    benchmark:
        'rsrc/output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.cache.rsrc',
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    run:
        import pandas as pd
        import numpy as np

        chrom_sizes = {}
        with open(input.fai, 'r') as faidx:
            for line in faidx:
                name, size = line.split()[:2]
                chrom_sizes[name] = int(size)

        with open(input.bg, 'r') as bedgraph:
            for line in bedgraph:
                chrom, start, end, cov = line.strip().split()
                last_chrom = chrom
                chrom_cov = np.zeros(chrom_sizes[chrom], dtype=np.int16)
                break

        with pd.HDFStore(output.hdf, mode='w', complib='blosc', complevel=9) as hdf:
            with open(input.bg, 'r') as bedgraph:
                for line in bedgraph:
                    chrom, start, end, cov = line.strip().split()
                    start = int(start) - 1
                    end = int(end)
                    if chrom != last_chrom:
                        out_key = f'{chrom}'
                        hdf.put(out_key, pd.Series(chrom_cov, dtype=np.int16), format='fixed')
                        last_chrom = chrom
                        chrom_cov = np.zeros(chrom_sizes[chrom], dtype=np.int16)
                    chrom_cov[start:end] = int(cov)
    # END OF RUN BLOCK


# Prepare meryl DBs for merqury-like QV estimation

rule meryl_query_only_kmer_db:
    """
    Create DB containing k-mers unique to the query sequences
    (the sequences for which the QV estimate should be computed)
    """
    input:
        query_db = 'output/kmer_db/{sample}_{readset1}.meryl',
        reference_db = 'output/kmer_db/{sample}_{readset2}.meryl'
    output:
        query_only = directory('output/kmer_db/{sample}_{readset1}_DIFF_{readset2}.meryl')
    benchmark:
        'rsrc/output/kmer_db/{sample}_{readset1}_DIFF_{readset2}.meryl.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * (attempt ** 4),
    shell:
        'meryl difference output {output.query_only} {input.query_db} {input.reference_db}'


def set_meryl_memory(wildcards, attempt):

    memory_mb = 8192 + 2048 * attempt
    if 'ONTUL' in wildcards.readset:
        if 'HIFIEC' in wildcards.readset and 'mbg' in wildcards.readset:
            pass
        else:
            # raw ONTUL reads
            memory_mb = 262144 + 65536 * attempt
    return memory_mb


rule meryl_generate_individual_kmer_stats:
    """
    -existence:
        Generate a tab-delimited line for each input sequence with the
        number of kmers in the sequence, in the database and common to both.
    """
    input:
        query_only = 'output/kmer_db/{sample}_{readset}_DIFF_{short_reads}.meryl',
        query_reads = select_winnowmap_reads
    output:
        table = 'output/kmer_stats/{sample}_{readset}_DIFF_{short_reads}.seqkm.tsv'
    benchmark:
        'rsrc/output/kmer_stats/{sample}_{readset}_DIFF_{short_reads}.seqkm.meryl.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: set_meryl_memory(wildcards, attempt),
        runtime_hrs = lambda wildcards, attempt: 8 * attempt
    shell:
        'meryl-lookup -existence -sequence {input.query_reads} -mers {input.query_only} > {output.table}'


def prob_base_correct(kmer_shared, kmer_total, kmer_size):
    return (kmer_shared / kmer_total) ** (1/kmer_size)


def base_error_rate(kmer_assembly_only, kmer_total, kmer_size):
    return 1 - (1 - kmer_assembly_only / kmer_total) ** (1/kmer_size)


def qv_estimate(error_rate):
    return -10 * math.log10(error_rate)


rule compute_global_query_qv_estimate:
    """
    Formulas for QV estimation as stated in 

    Rhie, A., Walenz, B.P., Koren, S. et al.
    Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies.
    Genome Biol 21, 245 (2020). https://doi.org/10.1186/s13059-020-02134-9

    Methods section "Consensus quality (QV) estimation"
    """
    input:
        query_stats = 'output/kmer_stats/{sample}_{readset1}.meryl.stats.h5',
        query_only_stats = 'output/kmer_stats/{sample}_{readset1}_DIFF_{readset2}.meryl.stats.h5'
    output:
        'output/qv_estimate/{sample}_{readset1}_REF_{readset2}.qv.tsv'
    run:
        import pandas as pd
        import math

        with pd.HDFStore(input.query_stats, mode='r') as hdf:
            num_kmer_query_total = hdf['statistics']['present']
            kmer_size_query_total = hdf['statistics']['kmer_size']

        with pd.HDFStore(input.query_only_stats, mode='r') as hdf:
            num_kmer_query_only = hdf['statistics']['present']
            kmer_size_query_only = hdf['statistics']['kmer_size']
        if kmer_size_query_total != kmer_size_query_only:
            raise ValueError(f'k-mer sizes do not match: {kmer_size_query_total} vs {kmer_size_query_only}')

        error_rate = base_error_rate(num_kmer_query_only, num_kmer_query_total, kmer_size_query_total)
        qv_est = round(-10 * math.log10(error_rate), 1)

        with open(output[0], 'w') as table:
            _ = table.write(f'sample\t{wildcards.sample}\n')
            _ = table.write(f'query_sequences\t{wildcards.readset1}\n')
            _ = table.write(f'reference_sequences\t{wildcards.readset2}\n')
            _ = table.write(f'kmer_size\t{kmer_size_query_only}\n')
            _ = table.write(f'kmer_query_only\t{num_kmer_query_only}\n')
            _ = table.write(f'kmer_query_total\t{num_kmer_query_total}\n')
            _ = table.write(f'error_rate\t{error_rate}\n')
            _ = table.write(f'QV_estimate\t{qv_est}\n')
    # END OF RUN BLOCK


rule compute_local_query_qv_estimate:
    """
    Comment regarding QV computation:
    The input meryl DB used for the lookup/existence operation contains
    only kmers unique to the (long read) sequences (i.e., the kmers from the
    short read dataset were removed [op difference]). Hence, the QV computation
    below uses "kmers shared between read sequence and read DB" as sequence-unique
    kmers (= not supported by Illumina/short read kmers), divided by the number
    of total kmers in the read sequence.
    """
    input:
        tsv = 'output/kmer_stats/{sample}_{readset}_DIFF_{short_reads}.seqkm.tsv'
    output:
        hdf = 'output/qv_estimate/{sample}_{readset}_REF_{short_reads}.seq-qv.h5'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    params:
        kmer_size = 31
    run:
        import pandas as pd
        import numpy as np

        table_header = ['read_name', 'kmer_in_seq', 'kmer_in_db', 'kmer_shared']
        df = pd.read_csv(input.tsv, sep='\t', names=table_header, header=None, index_col=None)

        # uncertain: for a test case of HiFi reads, the number of shared (= unsupported)
        # kmers was almost always 0. This might be an artifact due to the length
        # (statistical: expected number of errors per read for error rates ~10e-4 ~ 0),
        # or indicates some bug/problem...
        # For now, deal with that situation explicitly:

        df['error_rate'] = 0
        df['qv_estimate'] = -1

        select_nz = np.array(df['kmer_shared'] > 0, dtype=np.bool)

        df.loc[select_nz, 'error_rate'] = 1 - (1 - df.loc[select_nz, 'kmer_shared'] / df.loc[select_nz, 'kmer_in_seq']) ** (1/params.kmer_size)
        df.loc[select_nz, 'qv_estimate'] = (-10 * np.log10(df.loc[select_nz, 'error_rate'])).round(1)

        with pd.HDFStore(output.hdf, 'w', complevel=6) as hdf:
            hdf.put('sequence_qv', df, format='fixed')
