localrules: run_all


rule run_all:
    input:
        'output/alignments/ont_to_assm_graph/NA24385_giab_ULfastqs_guppy324_MAP-TO_v0151_2hap.r_utg.gaf',
        'output/alignments/ont_to_assm_graph/NA24385_ONT_PAD64459_Guppy32_MAP-TO_v0151_2hap.r_utg.gaf',
        'output/alignments/ont_to_assm_graph/NA24385_giab_ULfastqs_guppy324_MAP-TO_v0152_patg.r_utg.gaf',
        'output/alignments/ont_to_assm_graph/NA24385_ONT_PAD64459_Guppy32_MAP-TO_v0152_patg.r_utg.gaf',
        'output/alignments/ont_to_mbg_graph/NA24385_giab_ULfastqs_guppy324_MAP-TO_HIFIec_k2001_w1000.mbg.gaf',
        'output/alignments/ont_to_mbg_graph/NA24385_ONT_PAD64459_Guppy32_MAP-TO_HIFIec_k2001_w1000.mbg.gaf',


rule ont_error_correction:
    """
    """
    input:
        graph = 'input/{graph_type}_graph/{sample}.{graph}.{tigs}.gfa,
        reads = 'input/ont/{sample}_{readset}.fa.gz',
    output:
        gaf = 'output/alignments/ont_to_{graph_type}_graph/{sample}_{readset}_MAP-TO_{graph}.{tigs}.gaf',
        ec_reads_clip = 'output/alignments/ont_to_{graph_type}_graph/{sample}_{readset}_MAP-TO_{graph}.{tigs}.clip-ec.fa.gz',
        ec_reads_raw = 'output/alignments/ont_to_{graph_type}_graph/{sample}_{readset}_MAP-TO_{graph}.{tigs}.raw-ec.fa.gz',
    log:
        'log/output/alignments/ont_to_{graph_type}_graph/{sample}_{readset}_MAP-TO_{graph}.{tigs}.ga.log'
    benchmark:
        'rsrc/output/alignments/ont_to_{graph_type}_graph/{sample}_{readset}_MAP-TO_{graph}.{tigs}.ga.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 65536 + 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    shell:
        'GraphAligner -g {input.graph} -f {input.reads} '
            '-x vg -t {threads} '
            '--min-alignment-score 100 --multimap-score-fraction 1 '
            '--corrected-clipped-out {output.ec_reads_clip} --corrected-out {output.ec_reads_raw} '
            '-a {output.gaf} &> {log}'


rule build_hifi_read_graph:
    """
    old call for MBG v1.03:
    MBG -i {input} -o {output} -t {threads} --blunt -k {wildcards.kmer} -w {wildcards.window} &> {log}

    MBG v1.04+ contains fix for overlap-longer-than-nodes error

    MBG v1.05+ has better mem management, should limit <100G
    """
    input:
        'input/ec_hifi/NA24385_hpg_pbsq2-ccs_1000.ec-reads.fasta.gz'
    output:
        'input/mbg_graph/NA24385.HIFIec_k2001_w1000.mbg.gfa'
    log:
        'log/input/mbg_graph/NA24385.HIFIec_k2001_w1000.mbg.log'
    benchmark:
        'rsrc/input/mbg_graph/NA24385.HIFIec_k2001_w1000.mbg.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 65536 + 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    shell:
        'MBG -i {input} -o {output} -t {threads} -k 2001 -w 1000 &> {log}'
