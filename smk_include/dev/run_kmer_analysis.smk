
include: 'prep_custom_references.smk'

KMER_CONFIG = {
    'kmer_size': 31,  # increasing that requires changing the Bifrost build command!
    'trim_min_qual': 20,
    'trim_min_length': 51,
    'genome_size': int(3.1e9),  # this is roughly the sequence length of the HGSVC2 reference
    'ref_assembly': 'GRCh38_HGSVC2_noalt',
    'ref_annotation': 'ENSEMBLv98_RegBuild',
    'kmer_ratio': 99  # values below ~95 lower chances of detecting hap-specific sequence
}


def find_sample_short_reads(sample):

    data_sources = config['sample_description_' + sample]['data_sources']
    short_reads = []
    for readset_record in data_sources:
        for readset_type, readset_desc in readset_record.items():
            if readset_type != 'short_reads':
                continue
            readset_sample, readset_name = readset_desc['readset'].split('_', 1)
            assert readset_sample == sample, 'Sample mismatch: {} / {}'.format(sample, readset_desc)
            short_reads.append(readset_name)
    if not short_reads:
        raise ValueError('No short read data available for sample {}'.format(sample))
    return short_reads


def determine_possible_computations(wildcards):
    """
    NA19239_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta
    """
    module_outputs = {
        'annotation_table': 'output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.{annotation}.{ratio}.tsv',
        'kmer_counts': 'output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.kmer-counts.tsv',
    }

    fix_wildcards = {
        'ratio': KMER_CONFIG['kmer_ratio'],
        'annotation': KMER_CONFIG['ref_annotation'],
        'known_ref': KMER_CONFIG['ref_assembly']
    }

    compute_results = set()

    search_path = os.path.join(os.getcwd(), 'output/evaluation/kmer_analysis/phased_assemblies')
    for ps_assm in os.listdir(search_path):
        if not ps_assm.endswith('.fasta'):
            continue
        assm_base, hap, polisher, ext = ps_assm.split('.')
        sample, assm_reads = assm_base.split('_', 1)
        tmp = dict(fix_wildcards)
        tmp['sample'] = sample
        tmp['polisher'] = polisher
        tmp['assembly'] = assm_reads
        short_reads = find_sample_short_reads(sample)
        for sr in short_reads:
            tmp['readset'] = sr
            for target in module_outputs.values():
                fmt_target = target.format(**tmp)
                compute_results.add(fmt_target)
    
    return sorted(compute_results)


localrules: master_kmer_analysis,
            write_bifrost_fofn


rule master_kmer_analysis:
    input:
        determine_possible_computations


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
       'log/output/check_files/src_build/install_bifrost.log'
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
         'git checkout ab43065337a4f7c8ff9af9f77e6b3f2fe282b09f && '
         'mkdir build && cd build && '
         'cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX .. && '
         'make && make install ; ) > {log} 2>&1'


rule install_source_venn_diagram:
    input:
        'output/check_files/src_build/install_bifrost.ok'
    output:
        touch('output/check_files/src_build/install_vennd.ok')
    log:
       'log/output/check_files/src_build/install_vennd.log'
    conda:
        '../../environment/conda/conda_compile.yml'
    params:
        repo_folder = 'output/repositories'
    shell:
         '( rm -rf {params.repo_folder}/theupsetkmer && '
         'mkdir -p {params.repo_folder} && '
         'cd {params.repo_folder} && '
         'git clone https://github.com/tobiasmarschall/theupsetkmer.git && '
         'cd theupsetkmer && '
         'git checkout 58308d195f335ce2bc84747a5fd44bf37cd68656 && '
         'chmod u+x compile.sh && '
         './compile.sh && '
         'cp venn_diagram $CONDA_PREFIX/bin ; ) > {log} 2>&1'


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
        quality_trim = KMER_CONFIG['trim_min_qual'],
        min_read_length = KMER_CONFIG['trim_min_length'],
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
        runtime_hrs = lambda wildcards, attempt: attempt * 23,
        mem_per_cpu_mb = lambda wildcards, attempt: int(24576 * attempt / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 24676 * attempt
    params:
        kmer_size = KMER_CONFIG['kmer_size'],
        alpha = lambda wildcards, input: compute_lighter_alpha(input.report1, input.report2, KMER_CONFIG['genome_size']),
        genomesize = KMER_CONFIG['genome_size'],
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
        seq_mito = 'references/assemblies/{known_ref}.chrM.fasta'
    output:
        read_fofn = 'output/evaluation/kmer_analysis/{known_ref}/{sample}/{readset}.{assembly}.{polisher}.reads.txt',
        assm_fofn = 'output/evaluation/kmer_analysis/{known_ref}/{sample}/{readset}.{assembly}.{polisher}.assm.txt',
    run:
        import os
        with open(output.read_fofn, 'w') as dump:
            _ = dump.write(os.path.abspath(input.mate1) + '\n')
            _ = dump.write(os.path.abspath(input.mate2) + '\n')

        with open(output.assm_fofn, 'w') as dump:
            _ = dump.write(os.path.abspath(input.hap1) + '\n')
            _ = dump.write(os.path.abspath(input.hap2) + '\n')
            _ = dump.write(os.path.abspath(input.reference) + '\n')
            _ = dump.write(os.path.abspath(input.seq_mito) + '\n')


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
        kmer_size = KMER_CONFIG['kmer_size'],
        out_prefix = lambda wildcards, output: output[0].rsplit('.', 1)[0]
    shell:
        'Bifrost build --input-seq-file {input.read_fofn} --input-ref-file {input.assm_fofn} '
        '--output-file {params.out_prefix} --threads {threads} --colors --kmer-length {params.kmer_size} '
        '--verbose &> {log}'


rule count_kmers_per_color:
    """
    Usage: ./venn_diagram <graph_file.gfa> <color_file.bfg_colors> <k> <threads> <output_file.txt>
    """
    input:
        setup_ok = 'output/check_files/src_build/install_vennd.ok',
        graph = 'output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.gfa',
        colors = 'output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.bfg_colors'
    output:
        'output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.kmer-counts.tsv',
    log:
       'log/output/evaluation/kmer_analysis/{known_ref}/{sample}.{readset}.{assembly}.{polisher}.vd-count.log',
    benchmark:
        os.path.join('run/output/evaluation/kmer_analysis/{known_ref}',
                     '{sample}.{readset}.{assembly}.{polisher}.vd-count' + '.t{}.rsrc'.format(config['num_cpu_low']))
    conda: '../../environment/conda/conda_compile.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 16768 + 16768 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: int((16768 + 16768 * attempt) / config['num_cpu_low']),
        runtime_hrs = lambda wildcards, attempt: 2 * attempt
    params:
        kmer_size = KMER_CONFIG['kmer_size']
    shell:
        'venn_diagram {input.graph} {input.colors} {params.kmer_size} {threads} {output} &> {log}'


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
        mem_total_mb = lambda wildcards, attempt: 16768 + 16768 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: int((16768 + 16768 * attempt) / config['num_cpu_high']),
        runtime_hrs = lambda wildcards, attempt: 2 * attempt
    params:
        out_prefix = lambda wildcards, output: output[0].rsplit('.', 1)[0],
        kmer_ratio = lambda wildcards: round(float(wildcards.ratio) / 100, 2)
    shell:
         'Bifrost query --input-graph-file {input.graph} --input-color-file {input.colors} '
         '--input-query-file {input.queries} --output-file {params.out_prefix} '
         '--threads {threads} --inexact --ratio-kmers {params.kmer_ratio}'
         '--verbose &> {log}'
