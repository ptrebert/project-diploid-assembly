"""
Follow HiFi/ONT hybrid approach developed by MR
https://github.com/maickrau/hybrid-assembly/blob/master/commands.sh
[lines 46-65]
"""


rule hybrid_ga_align_ont_to_string_graph:
    """
    This rule uses a source built from GA's MultiseedDP branch
    to enable command line parameter "discard-cigar" (less relevant)
    and the computation of MAPQ scores (relevant, used in downstream scripts)
    """
    input:
        container = ancient('graphaligner.MultiseedDP.sif'),
        graph = lambda wildcards: ASSEMBLED_SAMPLES[wildcards.sample][wildcards.tigs],
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
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 176128 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt,
    shell:
        'module load Singularity && singularity exec {input.container} '
        'GraphAligner -g {input.graph} -f {input.reads} '
            '-t {threads} '
            '--min-alignment-score 5000 --multimap-score-fraction 0.99 '
            '--precise-clipping 0.7 --hpc-collapse-reads '
            '--seeds-mxm-length 30 --seeds-mem-count 10000 '
            '-b 15 --discard-cigar '
            '--corrected-out {output.hybrid_reads} '
            '-a {output.gaf} &> {log}'
