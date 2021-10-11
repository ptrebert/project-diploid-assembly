
################################
# This module is deprecated and  
# dysfunctional.
# Rules will be distributed to
# more task-specific modules
################################

def compute_read_coverage_files():

    out_path = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.Q{minq}.L{minl}.WS{window_size}.h5'

    sample=['NA18989']
    readset=[
        'ONTUL_guppy-5.0.11-sup-prom',
        'HIFIEC_hifiasm-v0.15.4',
        'ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k1001-w500',
        'ONTUL_guppy-5.0.11-sup-prom_MAP-TO_HIFIEC.mbg-k3001-w2000',
    ]
    reference=['T2Tv11_38p13Y_chm13']
    ql_combs = [(0, 0), (60, 0), (60, 20000), (60, 100000)]
    window_sizes = [20000, 50000, 100000]

    read_cov_files = []

    for s in sample:
        for rs in readset:
            for ref in reference:
                for (q, l) in ql_combs:
                    for ws in window_sizes:
                        formatter = {
                            'sample': s,
                            'readset': rs,
                            'reference': ref,
                            'minq': q,
                            'minl': l,
                            'window_size': ws
                        }
                        tmp = out_path.format(**formatter)
                        read_cov_files.append(tmp)
    return read_cov_files


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
        cache_read_cov = compute_read_coverage_files(),
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



rule cache_read_coverage:
    """
    input bed:
    chr1    1       16514   m64039_210506_192609/37554193/ccs       60      +
    """
    input:
        bed = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.aln.bed',
        fai = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta.fai'),
    output:
        hdf = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.Q{minq}.L{minl}.cache.h5',
    benchmark:
        'rsrc/output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.Q{minq}.L{minl}.cache.rsrc',
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt
    run:
        import pandas as pd
        import numpy as np
        import collections as col

        min_mapq = int(wildcards.minq)
        min_length = int(wildcards.minl)

        chrom_sizes = {}
        with open(input.fai, 'r') as faidx:
            for line in faidx:
                name, size = line.split()[:2]
                chrom_sizes[name] = int(size)

        with open(input.bed, 'r') as bed:
            chrom = bed.readline().split()[0]
            assert '#' not in chrom
            last_chrom = chrom
            chrom_cov = np.zeros(chrom_sizes[chrom], dtype=np.int16)

        mapqs = col.Counter()
        lengths = []
        aln = col.Counter() 

        with pd.HDFStore(output.hdf, mode='w', complib='blosc', complevel=9) as hdf:
            with open(input.bed, 'r') as bed:
                for line in bed:
                    chrom, start, end, read, mapq, strand = line.strip().split()
                    start = int(start) - 1
                    end = int(end)
                    mapq = int(mapq)
                    length = end - start
                    if chrom != last_chrom:
                        out_key = f'{last_chrom}/cov'
                        hdf.put(out_key, pd.Series(chrom_cov, dtype=np.int16), format='fixed')
                        out_key = f'{last_chrom}/mapq'
                        hdf.put(out_key, pd.Series(mapqs, dtype=np.int32), format='fixed')
                        out_key = f'{last_chrom}/aln_length'
                        hdf.put(out_key, pd.Series(lengths, dtype=np.int32), format='fixed')
                        out_key = f'{last_chrom}/reads'
                        read_abundance = pd.DataFrame.from_records(
                            [(k, v) for k, v in aln.items()],
                            columns=['read', 'num_aln']
                        )
                        hdf.put(out_key, read_abundance, format='fixed')
                        # reset
                        last_chrom = chrom
                        chrom_cov = np.zeros(chrom_sizes[chrom], dtype=np.int16)
                        mapqs = col.Counter()
                        lengths = []
                        aln = col.Counter()
                    if mapq < min_mapq or length < min_length:
                        continue
                    chrom_cov[start:end] += 1
                    mapqs[mapq] += 1
                    lengths.append(length)
                    aln[read] += 1

            out_key = f'{chrom}/cov'
            hdf.put(out_key, pd.Series(chrom_cov, dtype=np.int16), format='fixed')
            out_key = f'{chrom}/mapq'
            hdf.put(out_key, pd.Series(mapqs, dtype=np.int32), format='fixed')
            out_key = f'{chrom}/aln_length'
            hdf.put(out_key, pd.Series(lengths, dtype=np.int32), format='fixed')
            out_key = f'{chrom}/reads'
            read_abundance = pd.DataFrame.from_records(
                [(k, v) for k, v in aln.items()],
                columns=['read', 'num_aln']
            )
            hdf.put(out_key, read_abundance, format='fixed')
    # END OF RUN BLOCK


rule bin_read_coverage:
    input:
        hdf = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.Q{minq}.L{minl}.cache.h5',
    output:
        hdf = 'output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.Q{minq}.L{minl}.WS{window_size}.h5',
    benchmark:
        'rsrc/output/alignments/reads_to_linear_ref/{sample}_{readset}_MAP-TO_{reference}.Q{minq}.L{minl}.WS{window_size}.rsrc',
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt
    run:
        import pandas as pd
        import numpy as np

        ws = int(wildcards.window_size)
        with pd.HDFStore(input.hdf, 'r') as hdf_in:
            for key in hdf_in.keys():
                if not key.endswith('/cov'):
                    continue
                chrom = key.strip('/').split('/')[0]
                chrom_cov = hdf_in[key]
                if chrom_cov.size < ws:
                    # NB: chrM is part of T2T reference
                    continue
                blunt_end = chrom_cov.size // ws * ws

                window_avg_cov = chrom_cov.values[:blunt_end].reshape((-1, ws)).mean(axis=1)
                window_med_cov = np.sort(chrom_cov.values[:blunt_end].reshape((-1, ws)), axis=1)[:, ws//2]

                df = pd.DataFrame(
                    [window_avg_cov, window_med_cov],
                    index=['avg_cov', 'med_cov']
                )
                df = df.transpose()
                df['avg_cov_pctrk'] = df['avg_cov'].rank(ascending=True, pct=True)
                df['med_cov_pctrk'] = df['avg_cov'].rank(ascending=True, pct=True)

                with pd.HDFStore(output.hdf, 'a', complevel=9, complib='blosc') as hdf_out:
                    hdf_out.put(f'{chrom}/window_cov', df, format='fixed')
    # END OF RUN BLOCK


# Prepare meryl DBs for merqury-like QV estimation
