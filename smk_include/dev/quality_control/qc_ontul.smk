
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

rule run_all:
    input:
        INPUT_READS,
        'output/cdbg/NA18989.gfa',


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
