
rule build_hifi_read_dbg:
    input:
        hifi_ec_reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.read_type]
    output:
        graph = 'output/mbg_hifi/{sample}_{read_type}_{readset}.MBG-k{kmer}-w{window}.gfa',
        paths = 'output/mbg_hifi/{sample}_{read_type}_{readset}.MBG-k{kmer}-w{window}.gaf'
    log:
        'log/output/mbg_hifi/{sample}_{read_type}_{readset}.MBG-k{kmer}-w{window}.MBG.log'
    benchmark:
        'rsrc/output/mbg_hifi/{sample}_{read_type}_{readset}.MBG-k{kmer}-w{window}.MBG.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 65536 + 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 8 * attempt
    params:
        validate = lambda wildcards, input: validate_readset(wildcards.readset, input.hifi_ec_reads)
    shell:
        'MBG -i {input.hifi_ec_reads} -t {threads} '
            '-k {wildcards.kmer} -w {wildcards.window} '
            '--error-masking collapse-msat --include-end-kmers '
            '--out {output.graph} --output-sequence-paths {output.paths} &> {log}'


rule ont_to_graph_alignment:
    """
    """
    input:
        sif = ancient('graphaligner.MultiseedClusters.sif'),
        graph = 'output/mbg_hifi/{sample}_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.gfa',
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample]['ONTUL']
    output:
        gaf = 'output/alignments/ont_to_mbg_graph/{sample}_ONTUL_{readset}_MAP-TO_{graph_reads}.{graph_readset}.MBG-k{kmer}-w{window}.gaf',
        ec_reads_clip = 'input/ONTEC/{sample}_ONTEC_{readset}_MAP-TO_{graph_reads}.{graph_readset}.MBG-k{kmer}-w{window}.fasta.gz',
        ec_reads_full = 'input/ONTHY/{sample}_ONTHY_{readset}_MAP-TO_{graph_reads}.{graph_readset}.MBG-k{kmer}-w{window}.fasta.gz'
    log:
        'log/output/alignments/ont_to_mbg_graph/{sample}_ONTUL_{readset}_MAP-TO_{graph_reads}.{graph_readset}.MBG-k{kmer}-w{window}.GA.log'
    benchmark:
        'rsrc/output/alignments/ont_to_mbg_graph/{sample}_ONTUL_{readset}_MAP-TO_{graph_reads}.{graph_readset}.MBG-k{kmer}-w{window}.GA.rsrc'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 90112 + 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: 72 * attempt
    params:
        preset = 'dbg'
    shell:
        'module load Singularity && singularity exec '
        '--bind /:/hilbert {input.sif} '
        'GraphAligner -g {input.graph} -f /hilbert/{input.reads} '
            '-x {params.preset} -t {threads} '
            '--min-alignment-score 5000 --multimap-score-fraction 1 '
            '--hpc-collapse-reads '
            '--corrected-clipped-out {output.ec_reads_clip} '
            '--corrected-out {output.ec_reads_full} '
            '-a {output.gaf} &> {log}'
