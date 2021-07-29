
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

rule run_all:
    input:
        readsets = INPUT_READS,
        kmer_dbs = expand('output/kmer_db/{readset}.total.count', readset=READSETS),
        ont_align = expand(
            'output/alignments/ont_to_mbg_hifi/{sample}_{readset}_MAP-TO_HIFIEC.mbg-k{kmer}-w{window}.gaf',
            zip,
            sample=['NA18989'] * 4,
            readset=['ONTUL_guppy-5.0.11-sup-prom', 'ONTUL_guppy-4.0.11-hac-prom'] * 2,
            kmer=[1001, 3001] * 2,
            window=[500, 2000] * 2
        ),
        'output/kmer_stats/NA18989_ONTUL_guppy-5.0.11-sup-prom.unsupported.txt',
        'output/kmer_stats/NA18989_ONTUL_guppy-4.0.11-hac-prom.unsupported.txt'

        #'output/cdbg/NA18989.gfa',


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
    assert len(readset) == 1
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
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 90112 + 16384 * attempt,
        runtime_hrs = lambda wildcards, attempt: 24 * attempt
    params:
        preset = 'dbg'
    shell:
        'module load Singularity && singularity exec {input.container} '
        'GraphAligner -g {input.graph} -f {input.reads} '
            '-x {params.preset} -t {threads} '
            '--min-alignment-score 10000 --multimap-score-fraction 1 '
            '--corrected-clipped-out {output.ec_reads_clip} '
            '-a {output.gaf} &> {log}'


rule dump_unsupported_kmers:
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
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
    shell:
        'meryl print [difference {input.ont_db} {input.hifiec} {input.hifiaf} {input.short}] > {output}'
