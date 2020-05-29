
include: '../constraints.smk'
include: '../aux_utilities.smk'
include: 'prep_custom_references.smk'

localrules: master_kmer_analysis,
            write_bifrost_fofn


rule master_kmer_analysis:
    input:
        []

rule create_conda_environment_compile:
    output:
        'output/check_files/environment/conda_compile.ok'
    log:
        'log/output/check_files/environment/conda_compile.log'
    params:
        script_exec = lambda wildcards: find_script_path('inspect_environment.py', 'utilities')
    conda:
        '../../environment/conda/conda_compile.yml'
    shell:
        '{params.script_exec} '
        '--export-conda-env --outfile {output} --logfile {log}'


rule install_source_bifrost:
    """
    Bioconda version of Bifrost is built with AVX2 enabled,
    which does not work on older machines. Additionally, recent
    bug fix for graph querying has not been part of any release:
    github.com/pmelsted/bifrost/issues/20
    github.com/pmelsted/bifrost/issues/21
    """
    input:
        'output/check_files/environment/conda_compile.ok'
    output:
        touch('output/check_files/src_build/install_bifrost.ok')
    log:
       'output/check_files/src_build/install_bifrost.ok'
    conda:
        '../../environment/conda/conda_compile.yml'
    params:
        repo_folder = 'output/repositories'
    shell:
         '( rm -rf {params.repo_folder}/bifrost && '
         'mkdir -p {params.repo_folder} && '
         'cd {params.repo_folder} && '
         'git clone https://github.com/pmelsted/bifrost.git && '
         'cd bifrost && '
         'git checkout f1e48443f3576429590b65d24b998c88a52fb6d4 && '
         'mkdir build && cd build && '
         'cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX .. && '
         'make && make install ; ) > {log} 2>&1'


rule short_read_quality_trimming:
    """
    Note that due to the rather "hidden" way trim-galore/cutadapt
    are handling multithreading, the number of threads "4" is
    internally translated to something like 15 or 16 (reading, writing etc.)
    """
    input:
        mate1 = 'input/fastq/{readset}_1.fastq.gz',
        mate2 = 'input/fastq/{readset}_2.fastq.gz'
    output:
        mate1 = 'input/fastq/{readset}/trimmed/{readset}_1_val_1.fq.gz',
        report1 = 'input/fastq/{readset}/trimmed/{readset}_1.fastq.gz_trimming_report.txt',
        mate2 = 'input/fastq/{readset}/trimmed/{readset}_2_val_2.fq.gz',
        report2 = 'input/fastq/{readset}/trimmed/{readset}_2.fastq.gz_trimming_report.txt'
    log:
        'log/input/fastq/{readset}.trimming.log'
    benchmark:
        'run/input/fastq/{readset}.trimming.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: 16
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * 12,
        mem_per_cpu_mb = lambda wildcards, attempt: int(12288 * attempt / 16),
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
    params:
        quality_trim = 20,
        min_read_length = 51,
        outdir = lambda wildcards, input: os.path.join('input', 'fastq', wildcards.readset, 'trimmed')
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
        mate1 = 'input/fastq/{readset}/trimmed/{readset}_1_val_1.fq.gz',
        report1 = 'input/fastq/{readset}/trimmed/{readset}_1.fastq.gz_trimming_report.txt',
        mate2 = 'input/fastq/{readset}/trimmed/{readset}_2_val_2.fq.gz',
        report2 = 'input/fastq/{readset}/trimmed/{readset}_2.fastq.gz_trimming_report.txt'
    output:
        mate1 = 'input/fastq/{readset}/corrected/{readset}_1_val_1.cor.fq.gz',
        mate2 = 'input/fastq/{readset}/corrected/{readset}_2_val_2.cor.fq.gz',
    log:
        'log/input/fastq/{readset}.corr.log'
    benchmark:
        'run/input/fastq/{readset}.corr.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * 12,
        mem_per_cpu_mb = lambda wildcards, attempt: int(24576 * attempt / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 24676 * attempt
    params:
        kmer_size = 31,
        alpha = lambda wildcards, input: compute_lighter_alpha(input.report1, input.report2, 3272100381),
        genomesize = 3272100381,
        outdir = lambda wildcards, output: os.path.dirname(output.mate1)
    shell:
        'lighter -r {input.mate1} -r {input.mate2} '
        '-k {params.kmer_size} {params.genomesize} {params.alpha} '
        '-od {params.outdir} -t {threads} -zlib 6 &> {log}'


rule write_bifrost_fofn:
    input:
        mate1 = 'input/fastq/{sample}_{readset}/corrected/{sample}_{readset}_1_val_1.cor.fq.gz',
        mate2 = 'input/fastq/{sample}_{readset}/corrected/{sample}_{readset}_2_val_2.cor.fq.gz',
        hap1 = 'output/evaluation/kmer_analysis/phased_assemblies/{sample}_{assembly}.h1-un.{polisher}.fasta',
        hap2 = 'output/evaluation/kmer_analysis/phased_assemblies/{sample}_{assembly}.h2-un.{polisher}.fasta',
        reference = 'references/assemblies/{known_ref}.no-mito.fasta',
    output:
        read_fofn = 'output/evaluation/kmer_analysis/{known_ref}/{sample}/{readset}.{assembly}.{polisher}.reads.txt',
        assm_fofn = 'output/evaluation/kmer_analysis/{known_ref}/{sample}/{readset}.{assembly}.{polisher}.assm.txt',
    run:
        import os
        with open(output.read_fofn, 'w') as dump:
            _ = dump.write(os.path.abspath(input.mate1) + '\n')
            _ = dump.write(os.path.abspath(input.mate2))

        with open(output.assm_fofn, 'w') as dump:
            _ = dump.write(os.path.abspath(input.hap1) + '\n')
            _ = dump.write(os.path.abspath(input.hap2) + '\n')
            _ = dump.write(os.path.abspath(input.reference))


rule build_bifrost_colored_dbg:
    input:
        setup_ok = 'output/check_files/src_build/install_bifrost.ok',
        read_fofn = 'output/evaluation/kmer_analysis/{known_ref}/{sample}/{readset}.{assembly}.{polisher}.reads.txt',
        assm_fofn = 'output/evaluation/kmer_analysis/{known_ref}/{sample}/{readset}.{assembly}.{polisher}.assm.txt',
    output:
        'output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.gfa',
        'output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.bfg_colors'
    log:
       'log/output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.build.log',
    benchmark:
        os.path.join('run/output/evaluation/kmer_analysis/{known_ref}',
                     '{sample}.{readset}.{assembly}.{polisher}.build' + '.t{}.rsrc'.format(config['num_cpu_high']))
    conda: '../../environment/conda/conda_compile.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: int((32768 + 32768 * attempt) / config['num_cpu_high']),
        runtime_hrs = lambda wildcards, attempt: 16 * attempt
    params:
        kmer_size = 31,
        out_prefix = lambda wildcards, output: output[0].rsplit('.', 1)[0]
    shell:
        'Bifrost build --input-seq-file {input.read_fofn} --input-ref-file {input.assm_fofn} '
        '--output-file {params.out_prefix} --threads {threads} --colors --kmer-length {params.kmer_size} '
        '--verbose &> {log}'


rule query_bifrost_colored_dbg:
    input:
        setup_ok = 'output/check_files/src_build/install_bifrost.ok',
        graph = 'output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.gfa',
        colors = 'output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.bfg_colors',
        queries = 'references/annotation/{known_ref}-{annotation}.fasta'
    output:
        'output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.{annotation}.{ratio}.tsv',
    log:
        'log/output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.{annotation}.{ratio}.query.log',
    benchmark:
        os.path.join('run/output/evaluation/kmer_analysis/{known_ref}',
                     '{sample}.{readset}.{assembly}.{polisher}.{annotation}.{ratio}.query' + '.t{}.rsrc'.format(config['num_cpu_high']))
    conda: '../../environment/conda/conda_compile.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: int((32768 + 32768 * attempt) / config['num_cpu_high']),
        runtime_hrs = lambda wildcards, attempt: 12 * attempt
    params:
        out_prefix = lambda wildcards, output: output[0].rsplit('.', 1)[0],
        kmer_ratio = lambda wildcards: round(float(wildcards.ratio) / 100, 2)
    shell:
         'Bifrost query --input-graph-file {input.graph} --input-color-file {input.colors} '
         '--input-query-file {input.queries} --output-file {params.out_prefix} '
         '--threads {threads} --inexact --ratio-kmers {params.kmer_ratio}'
         '--verbose &> {log}'
