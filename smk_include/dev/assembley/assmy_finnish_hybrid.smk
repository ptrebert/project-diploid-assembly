"""
Follow HiFi/ONT hybrid approach developed by MR
https://github.com/maickrau/hybrid-assembly/blob/master/commands.sh
[lines 46-65]
"""


rule hybrid_ga_align_ont_to_string_graph:
    """
    This rule uses a source built from GA's MultiseedClusters branch
    to enable command line parameter "discard-cigar" (less relevant)
    and the computation of MAPQ scores (relevant, used in downstream scripts)

    NB: since the hifiasm graph used as input here is not HPC (but the
    string graph of the chrY T2T was), the alignment below omits the
    parameter "--hpc-collapse-reads" (cf. command list above)
    """
    input:
        container = ancient('graphaligner.MultiseedClusters.sif'),
        graph = 'output/clean_graphs/{sample_info}_{sample}.{tigs}.gfa',
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.ont_type]
    output:
        gaf = 'output/hybrid/ont_to_graph/{sample_info}_{sample}.{ont_type}.{tigs}.gaf',
        hybrid_reads = 'output/hybrid/ont_to_graph/{sample_info}_{sample}.{ont_type}.{tigs}.ONTHY.fasta',
    log:
        'log/output/hybrid/ont_to_graph/{sample_info}_{sample}.{ont_type}.{tigs}.ga.log',
    benchmark:
        'rsrc/output/hybrid/ont_to_graph/{sample_info}_{sample}.{ont_type}.{tigs}.ga.rsrc',
    wildcard_constraints:
        sample = CONSTRAINT_SAMPLES
#    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_max']
    resources:
        mem_total_mb = lambda wildcards, attempt: 1048576 * attempt,
        runtime_hrs = lambda wildcards, attempt: 167,
    shell:
        'module load Singularity && singularity exec '
        '--bind /:/hilbert {input.container} '
        'GraphAligner -g {input.graph} -f /hilbert/{input.reads} '
            '-t {threads} '
            '--min-alignment-score 5000 --multimap-score-fraction 0.99 '
            '--precise-clipping 0.7 '
            '--seeds-mxm-length 30 --seeds-mem-count 10000 '
            '-b 15 --discard-cigar '
            '--corrected-out {output.hybrid_reads} '
            '-a {output.gaf} &> {log}'


AFR_SAMPLE_1 = 'NA19317'
AFR_SAMPLE_2 = 'NA19347'

rule assemble_afr_mix:
    input:
        s1_reads = lambda wildcards: SAMPLE_INFOS[AFR_SAMPLE_1]['HIFIAF'],
        s2_reads = lambda wildcards: SAMPLE_INFOS[AFR_SAMPLE_2]['HIFIAF']
    output:
        primary_contigs = 'output/hybrid/afr_mix/hifiasm/{afr_sample1}-U-{afr_sample2}.p_ctg.gfa',
    log:
        hifiasm = 'log/output/hybrid/afr_mix/hifiasm/{afr_sample1}-U-{afr_sample2}.hifiasm.log',
    benchmark:
        'rsrc/output/hybrid/afr_mix/hifiasm/{afr_sample1}-U-{afr_sample2}.hifiasm.rsrc',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((180224 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 180224 * attempt,
        runtime_hrs = lambda wildcards, attempt: 36 * attempt
    params:
        prefix = lambda wildcards, output: output.primary_contigs.rsplit('.', 2)[0],
    shell:
        'hifiasm -o {params.prefix} -t {threads} --write-ec --write-paf --primary {input.s1_reads} {input.s2_reads} &> {log.hifiasm}'
