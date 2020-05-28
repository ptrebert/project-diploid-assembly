
include: '../aux_utilities.smk'

localrules: master_kmer_analysis


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
         'rm -rf {params.repo_folder}/bifrost && '
         'mkdir -p {params.repo_folder} && '
         'cd {params.repo_folder} && '
         'git clone https://github.com/pmelsted/bifrost.git && '
         'cd bifrost && '
         'git checkout f1e48443f3576429590b65d24b998c88a52fb6d4 && '
         'mkdir build && cd build && '
         'cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX .. && '
         'make && make install &> {log}'


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
        mem_per_cpu_mb = lambda wildcards, attempt: 12288 * attempt / 16,
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
    params:
        quality_trim = 20,
        min_read_length = 51,
        outdir = lambda wildcards, input: os.path.join('input', 'fastq', wildcards.readset, 'trimmed')
    shell:
        'trim_galore --quality {params.quality_trim} --length {params.min_read_length} '
        '--trim-n --output_dir {params.outdir} --cores 4 --paired '
        '{input.mate1} {input.mate2} &> {log}'


rule short_read_error_correction:
    input:
        mate1 = 'input/fastq/{readset}/trimmed/{readset}_1_val_1.fq.gz',
        mate2 = 'input/fastq/{readset}/trimmed/{readset}_2_val_2.fq.gz',
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
        mem_per_cpu_mb = lambda wildcards, attempt: 24576 * attempt / config['num_cpu_high'],
        mem_total_mb = lambda wildcards, attempt: 24676 * attempt
    params:
        kmer_size = 31,
        alpha = 0.1,
        outdir = lambda wildcards, output: os.path.dirname(output.mate1)
    shell:
        'lighter -r {input.mate1} -r {input.mate2} '
        '-k {params.kmer_size} {params.genomesize} {params.alpha} '
        '-od {params.outdir} -t {threads} -zlib 4 &> {log}'