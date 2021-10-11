import pathlib as pl

def set_errormasking(wildcards):

    if 'HIFIAF' in wildcards.read_type:
        errm = 'collapse-msat'
    elif 'HIFIEC' in wildcards.read_type:
        errm = 'no'
    else:
        raise ValueError(str(wildcards))
    return errm


def get_read_path(sample, read_type, prefix=''):

    file_path = str(SAMPLE_INFOS[sample][read_type])
    if prefix:
        file_path = str(pl.Path(prefix, file_path.strip('/')))
        assert 'gpfs' in file_path, f'Not an abs path: {file_path}'
    return file_path


rule build_hifi_read_dbg:
    """
    sif = ancient('mbg.master.sif'),
    """
    input:
        sif = ancient('mbg.UnitigResolve.sif'),
        reads = lambda wildcards: get_read_path(wildcards.sample, wildcards.read_type)
    output:
        graph = 'output/mbg_hifi/{sample}_{read_type}_{readset}.MBG-k{kmer}-w{window}.gfa',
        paths = 'output/mbg_hifi/{sample}_{read_type}_{readset}.MBG-k{kmer}-w{window}.gaf'
    log:
        'log/output/mbg_hifi/{sample}_{read_type}_{readset}.MBG-k{kmer}-w{window}.MBG.log'
    benchmark:
        'rsrc/output/mbg_hifi/{sample}_{read_type}_{readset}.MBG-k{kmer}-w{window}.MBG.rsrc'
    wildcard_constraints:
        read_type = '(HIFIEC|HIFIAF)'
#    conda:
#        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 65536 + 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    params:
        masking = set_errormasking,
        input_path = lambda wildcards: get_read_path(wildcards.sample, wildcards.read_type, '/hilbert')
    shell:
        'module load Singularity && singularity exec '
        '--bind /:/hilbert {input.sif} '
        'MBG -i {params.input_path} -t {threads} '
            '-k {wildcards.kmer} -w {wildcards.window} '
            '--error-masking {params.masking} --include-end-kmers '
            '--out {output.graph} --output-sequence-paths {output.paths} &> {log}'


def set_read_hpc(wildcards):

    if 'HIFIAF' in wildcards.graph_reads:
        hpc = '--hpc-collapse-reads'
    elif 'HIFIEC' in wildcards.graph_reads:
        hpc = ''
    else:
        raise ValueError(str(wildcards))
    return hpc


rule ont_to_graph_alignment:
    """
    """
    input:
        sif = ancient('graphaligner.MultiseedClusters.sif'),
        graph = 'output/mbg_hifi/{sample}_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.gfa',
        reads = 'input/ONTUL/{sample}_ONTUL_guppy-5.0.11-sup-prom.fasta.gz'
    output:
        gaf = 'output/alignments/ont_to_mbg_graph/{sample}_ONTUL_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.gaf',
        ec_reads_clip = 'output/alignments/ont_to_mbg_graph/{sample}_ONTEC_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.fasta.gz',
        ec_reads_full = 'output/alignments/ont_to_mbg_graph/{sample}_ONTHY_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.fasta.gz'
    log:
        'log/output/alignments/ont_to_mbg_graph/{sample}_ONTUL_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.GA.log'
    benchmark:
        'rsrc/output/alignments/ont_to_mbg_graph/{sample}_ONTUL_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.GA.rsrc'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 90112 + 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 72 * attempt
    params:
        preset = 'dbg',
        hpc = set_read_hpc,
    shell:
        'module load Singularity && singularity exec '
        '--bind /:/hilbert {input.sif} '
        'GraphAligner -g {input.graph} -f {input.reads} '
            '-x {params.preset} -t {threads} '
            '--min-alignment-score 5000 --multimap-score-fraction 1 '
            ' {params.hpc} '
            '--corrected-clipped-out {output.ec_reads_clip} '
            '--corrected-out {output.ec_reads_full} '
            '-a {output.gaf} &> {log}'


rule compute_ont_corrected_stats:
    input:
        'output/alignments/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.fasta.gz',
    output:
        'output/alignments/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.stats.tsv.gz',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 4
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: 12 * attempt
    shell:
        'pigz -p 2 -d -c {input} | seqtk comp | pigz -p 2 --best > {output} '


localrules: deduplicate_ont_corrected_reads

rule deduplicate_ont_corrected_reads:
    input:
        reads = expand(
            'output/alignments/ont_to_mbg_graph/{{sample}}_{{read_type}}_{{readset}}_MAP-TO_{{graph_reads}}_{{graph_readset}}.MBG-k{kmer}-w{window}.stats.tsv.gz',
            zip,
            kmer=config['mbg_kmers'],
            window=config['mbg_windows']
        )
    output:
        'output/alignments/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.foobar'
    run:
        raise