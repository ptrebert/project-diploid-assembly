
rule qc_merge_ontul:
    input:
        fastq = lambda wildcards: SAMPLE_INFOS[wildcards.sample]['ONTUL_RAW']
    output:
        'input/ONTUL/{sample}_ONTUL_' + f'{RS_ONTUL}.fasta.gz'
    benchmark:
        'rsrc/input/ONTUL/{sample}_ONTUL_' + f'{RS_ONTUL}.merge.rsrc'
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    shell:
        'pigz -d -c {input} | seqtk seq -A -C | pigz -p {threads} --best > {output}'


rule qc_local_seq_stats:
    input:
        'input/{read_type}/{sample}_{read_type}_{readset}.fasta.gz'
    output:
        'input/{read_type}/{sample}_{read_type}_{readset}.stats.tsv.gz'
    conda: '../../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        read_type = '(ONTUL|ONTEC|ONTHY)'
    threads: 2
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    shell:
        'pigz -d -c {input} | seqtk comp | pigz --best > {output}'


rule qc_remote_seq_stats:
    input:
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.read_type]
    output:
        'input/{read_type}/{sample}_{read_type}_{readset}.stats.tsv.gz'
    conda: '../../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        read_type = '(HIFIEC|HIFIAF|SHORT)'
    threads: 2
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    shell:
        'pigz -d -c {input} | seqtk comp | pigz --best > {output}'


rule qc_short_read_quality_trimming:
    """
    Note that due to the rather "hidden" way trim-galore/cutadapt
    are handling multithreading, the number of threads "4" is
    internally translated to something like 15 or 16 (reading, writing etc.)
    """
    input:
        mate1 = lambda wildcards: SAMPLE_INFOS[wildcards.sample]['SHORT_1'],
        mate2 = lambda wildcards: SAMPLE_INFOS[wildcards.sample]['SHORT_2'],
    output:
        mate1 = 'input/SHORT/{sample}/trimmed/{sample}_{readset}_1_val_1.fq.gz',
        report1 = 'input/SHORT/{sample}/trimmed/{sample}_{readset}_1.fastq.gz_trimming_report.txt',
        mate2 = 'input/SHORT/{sample}/trimmed/{sample}_{readset}_2_val_2.fq.gz',
        report2 = 'input/SHORT/{sample}/trimmed/{sample}_{readset}_2.fastq.gz_trimming_report.txt'
    log:
        'log/input/SHORT/{sample}_{readset}.trimming.log'
    benchmark:
        'run/input/SHORT/{sample}_{readset}.trimming.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 16
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * 12,
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
    params:
        quality_trim = 20,
        min_read_length = 51,
        outdir = lambda wildcards, input: os.path.join('input', 'SHORT', wildcards.sample, 'trimmed')
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


rule qc_short_read_error_correction:
    input:
        mate1 = 'input/SHORT/{sample}/trimmed/{sample}_{readset}_1_val_1.fq.gz',
        report1 = 'input/SHORT/{sample}/trimmed/{sample}_{readset}_1.fastq.gz_trimming_report.txt',
        mate2 = 'input/SHORT/{sample}/trimmed/{sample}_{readset}_2_val_2.fq.gz',
        report2 = 'input/SHORT/{sample}/trimmed/{sample}_{readset}_2.fastq.gz_trimming_report.txt'
    output:
        mate1 = 'input/SHORT/{sample}/corrected/{sample}_{readset}_1_val_1.cor.fq.gz',
        mate2 = 'input/SHORT/{sample}/corrected/{sample}_{readset}_2_val_2.cor.fq.gz',
    log:
        'log/input/SHORT/{sample}_{readset}.corr.log'
    benchmark:
        'run/input/SHORT/{sample}_{readset}.corr.rsrc'
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


rule qc_merge_shortec_reads:
    """
    This rule is more convenience s.t. all downstream rules can stay as-is
    """
    input:
        mate1 = 'input/SHORT/{sample}/corrected/{sample}_{readset}_1_val_1.cor.fq.gz',
        mate2 = 'input/SHORT/{sample}/corrected/{sample}_{readset}_2_val_2.cor.fq.gz',
    output:
        'input/SHORT/{sample}_SHORT_{readset}.fasta.gz'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: 24 * attempt
    shell:
        'pigz -d -c {input} | seqtk seq -A -C | pigz -p {threads} --best > {output}'
