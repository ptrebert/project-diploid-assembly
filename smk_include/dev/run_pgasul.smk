localrules: run_all


def find_script_path(script_name, subfolder=''):
    """
    Find full path to script to be executed. Function exists
    to avoid config parameter "script_dir"

    :param script_name:
    :param subfolder:
    :return:
    """
    import os

    current_root = workflow.basedir
    last_root = ''

    script_path = None

    for _ in range(workflow.basedir.count('/')):
        if last_root.endswith('project-diploid-assembly'):
            raise RuntimeError('Leaving project directory tree (next: {}). '
                               'Cannot find script {} (subfolder: {}).'.format(current_root, script_name, subfolder))
        check_path = os.path.join(current_root, 'scripts', subfolder).rstrip('/')  # if subfolder is empty string
        if os.path.isdir(check_path):
            check_script = os.path.join(check_path, script_name)
            if os.path.isfile(check_script):
                script_path = check_script
                break
        last_root = current_root
        current_root = os.path.split(current_root)[0]

    if script_path is None:
        raise RuntimeError('Could not find script {} (subfolder {}). '
                           'Started at path: {}'.format(script_name, subfolder, workflow.basedir))
    return script_path


rule run_all:
    input:
        'output/alignments/ont_to_assm_graph/NA24385_giab_ULfastqs_guppy324_MAP-TO_v0152_patg.r_utg.gaf',
        'output/alignments/ont_to_assm_graph/NA24385_ONT_PAD64459_Guppy32_MAP-TO_v0152_patg.r_utg.gaf',
        'output/alignments/ont_to_mbg_graph/NA24385_giab_ULfastqs_guppy324_MAP-TO_HIFIec_k2001_w1000.mbg.gaf',
        'output/alignments/ont_to_mbg_graph/NA24385_ONT_PAD64459_Guppy32_MAP-TO_HIFIec_k2001_w1000.mbg.gaf',
        #'output/alignments/ont_to_assm_graph/NA24385_giab_ULfastqs_guppy324_MAP-TO_v0151_2hap.r_utg.gaf',
        #'output/alignments/ont_to_assm_graph/NA24385_ONT_PAD64459_Guppy32_MAP-TO_v0151_2hap.r_utg.gaf',


wildcard_constraints:
    sample = 'NA24385'


rule check_overlong_edges:
    input:
        gfa = 'input/{graph_type}_graph/{sample}.{graph}.{tigs}.gfa'
    output:
        discard = 'input/{graph_type}_graph/{sample}.{graph}.{tigs}.discard.links'
    conda: '../../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('gfa_check_ovl.py'),
        stats_out = lambda wildcards, output: output.discard.replace('.discard.links', '.stats')
    resources:
        mem_total_mb = lambda wildcards, attempt: 512 * attempt,
        runtime_hrs = 0,
        runtime_min = lambda wildcards, attempt: 10 * attempt,
    shell:
        '{params.script_exec} -g {input.gfa} > {params.stats_out}'


rule clean_input_gfa:
    input:
        gfa = 'input/{graph_type}_graph/{sample}.{graph}.{tigs}.gfa',
        discard = 'input/{graph_type}_graph/{sample}.{graph}.{tigs}.discard.links'
    output:
        gfa = 'input/{graph_type}_graph/{sample}.{graph}.{tigs}.clean.gfa'
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    run:
        import os
        import io
        import shutil

        skip_lines = []
        with open(input.discard, 'r') as table:
            for line in table:
                if line.startswith('#') or not line.strip():
                    continue
                ln = line.split('\t')[0]
                skip_lines.append(int(ln))
        
        if not skip_lines:
            shutil.copy(input.gfa, output.gfa)
        else:
            gfa_buffer = io.StringIO()
            with open(input.gfa, 'r') as gfa:
                for ln, line in enumerate(gfa, start=1):
                    if ln in skip_lines:
                        continue
                    gfa_buffer.write(line)
            with open(output.gfa, 'w') as gfa:
                _ = gfa.write(gfa_buffer.getvalue())
    # END OF RUN BLOCK


rule ont_to_graph_alignment:
    """
    """
    input:
        container = ancient('graphaligner.sif'),
        graph = 'input/{graph_type}_graph/{sample}.{graph}.{tigs}.clean.gfa',
        reads = 'input/ont/{sample}_{readset}.fa.gz',
    output:
        gaf = 'output/alignments/ont_to_{graph_type}_graph/{sample}_{readset}_MAP-TO_{graph}.{tigs}.gaf',
        ec_reads_clip = 'output/alignments/ont_to_{graph_type}_graph/{sample}_{readset}_MAP-TO_{graph}.{tigs}.clip-ec.fa.gz',
        ec_reads_raw = 'output/alignments/ont_to_{graph_type}_graph/{sample}_{readset}_MAP-TO_{graph}.{tigs}.raw-ec.fa.gz',
    log:
        'log/output/alignments/ont_to_{graph_type}_graph/{sample}_{readset}_MAP-TO_{graph}.{tigs}.ga.log'
    benchmark:
        'rsrc/output/alignments/ont_to_{graph_type}_graph/{sample}_{readset}_MAP-TO_{graph}.{tigs}.ga.rsrc'
#    conda:
#        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 98304 + 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 72 * attempt
    params:
        preset = lambda wildcards: 'dbg' if wildcards.tigs == 'mbg' else 'vg'
    shell:
        'module load Singularity && singularity exec {input.container} '
        'GraphAligner -g {input.graph} -f {input.reads} '
            '-x {params.preset} -t {threads} '
            '--min-alignment-score 500 --multimap-score-fraction 1 '
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
