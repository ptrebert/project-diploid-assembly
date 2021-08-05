import pathlib

localrules: run_all, collect_aligned_strandseq, write_haploclust_config_json


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


def collect_sseq_libs_per_sample():

    import os

    lut = dict()

    for root, dirs, files in os.walk('input/sseq/'):
        if not files:
            continue
        libs = sorted([f.rsplit('_', 1)[-2] for f in files if f.endswith('_1.fastq.gz')])
        assert len(libs) > 1
        sample = os.path.split(root)[-1]
        assert all(sample in lib for lib in libs)
        lut[sample] = libs
    return lut


SSEQ_SAMPLE_TO_LIBS = collect_sseq_libs_per_sample()


rule run_all:
    input:
        'output/alignments/ont_to_assm_graph/NA24385_giab_ULfastqs_guppy324_MAP-TO_v0152_patg.r_utg.gaf',
        'output/alignments/ont_to_assm_graph/NA24385_ONT_PAD64459_Guppy32_MAP-TO_v0152_patg.r_utg.gaf',
        'output/alignments/ont_to_mbg_graph/NA24385_giab_ULfastqs_guppy324_MAP-TO_HIFIec_k2001_w1000.mbg.gaf',
        'output/alignments/ont_to_mbg_graph/NA24385_ONT_PAD64459_Guppy32_MAP-TO_HIFIec_k2001_w1000.mbg.gaf',
        'output/alignments/sseq_to_assm_graph/NA24385_MAP-TO_v0152_patg.r_utg.fofn',
        'output/alignments/sseq_to_mbg_graph/NA24385_MAP-TO_HIFIec_k2001_w1000.mbg.fofn',
        'output/haploclust/NA24385/sseq_to_assm_graph/NA24385_MAP-TO_v0152_patg.r_utg/haploclust.cfg.json',
        'output/haploclust/NA24385/sseq_to_mbg_graph/NA24385_MAP-TO_HIFIec_k2001_w1000.mbg/haploclust.cfg.json'


wildcard_constraints:
    sample = 'NA24385',
    tigs = '(r_utg|mbg)'


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
        mem_total_mb = lambda wildcards, attempt: 90112 + 16384 * attempt,
        runtime_hrs = lambda wildcards, attempt: 24 * attempt if wildcards.tigs == 'mbg' else 167
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


def extract_read_length(library_name):

    components = library_name.split('_')
    library_info = components[2].split('-')[1].strip('pe')
    read_length = int(library_info)
    return read_length


rule merge_strandseq_reads:
    input:
        mate1 = 'input/sseq/{sample}/{library_id}_1.fastq.gz',
        mate2 = 'input/sseq/{sample}/{library_id}_2.fastq.gz',
    output:
        merged = 'output/sseq_merged/temp/{sample}/{library_id}.assembled.fastq',
        discarded = 'output/sseq_merged/temp/{sample}/{library_id}.discarded.fastq',
        skip_frw = 'output/sseq_merged/temp/{sample}/{library_id}.unassembled.forward.fastq',
        skip_rev = 'output/sseq_merged/temp/{sample}/{library_id}.unassembled.reverse.fastq',
    log: 'log/output/sseq_merged/temp/{sample}/{library_id}.pear.log',
    benchmark: 'rsrc/output/sseq_merged/temp/{sample}/{library_id}.pear.rsrc',
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_pear_mb = lambda wildcards, attempt: 512 * attempt
    params:
        min_trim_length = lambda wildcards: int(round(extract_read_length(wildcards.library_id) * 0.8, 0)),
        output_prefix = lambda wildcards, output: output.merged.rsplit('.', 2)[0] 
    shell:
        'pear -f {input.mate1} -r {input.mate2} -t {params.min_trim_length} '
            '--memory {resources.mem_pear_mb}M -o {params.output_prefix} &> {log}'


rule concat_merged_and_unassembled_strandseq_reads:
    input:
        merged = 'output/sseq_merged/temp/{sample}/{library_id}.assembled.fastq',
        unassembled = 'output/sseq_merged/temp/{sample}/{library_id}.unassembled.forward.fastq',
    output:
        fa = 'output/sseq_merged/{sample}/{library_id}.merged.fasta.gz'
    benchmark: 'rsrc/output/sseq_merged/{sample}/{library_id}.merged.rsrc',
    threads: 3
    conda:
        '../../environment/conda/conda_biotools.yml'
    params:
        script_exec = lambda wildcards: find_script_path('fasta_add_libname.py'),
        compress_threads = 2
    shell:
        'seqtk seq -A -C {input.merged} {input.unassembled} | '
        '{params.script_exec} --lib-name {wildcards.library_id} | '
        'pigz -p {params.compress_threads} --best > {output.fa}'


rule dump_graph_to_fasta:
    input:
        graph = 'input/{graph_type}_graph/{sample}.{graph}.{tigs}.clean.gfa',
    output:
        fasta = 'input/{graph_type}_graph/{sample}.{graph}.{tigs}.clean.fasta',
        stats = 'input/{graph_type}_graph/{sample}.{graph}.{tigs}.clean.stats.tsv',
    threads: config['num_cpu_low']
    conda:
        '../../environment/conda/conda_pyscript.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 + 4096 * attempt,
    params:
        script_exec = lambda wildcards: find_script_path('gfa_to_fasta.py'),
    shell:
        '{params.script_exec} --graph {input.graph} --n-cpus {threads} '
            '--out-fasta {output.fasta} --out-stats {output.stats} --out-map /dev/null '
            '--break-lines 0 '


rule bwa_index_fasta:
    input:
        fasta = 'input/{graph_type}_graph/{sample}.{graph}.{tigs}.clean.fasta',
    output:
        multiext(
            'input/{graph_type}_graph/{sample}.{graph}.{tigs}.clean/bwa_index',
            '.amb', '.ann', '.bwt', '.pac', '.sa'
        )
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 + 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    params:
        prefix = lambda wildcards, output: output[0].rsplit('.', 1)[0]
    shell:
        'bwa index -p {params.prefix} {input.fasta}'


rule align_strandseq_to_tigs:
    input:
        index = 'input/{graph_type}_graph/{sample}.{graph}.{tigs}.clean/bwa_index.bwt',
        reads = 'output/sseq_merged/{sample}/{library_id}.merged.fasta.gz'
    output:
        bam = 'output/alignments/sseq_to_{graph_type}_graph/{sample}/{library_id}_MAP-TO_{graph}.{tigs}.psort.mdup.bam',
        bai = 'output/alignments/sseq_to_{graph_type}_graph/{sample}/{library_id}_MAP-TO_{graph}.{tigs}.psort.mdup.bam.bai',
        tmp_bam = temp('output/alignments/temp/sseq_to_{graph_type}_graph/{sample}/{library_id}_MAP-TO_{graph}.{tigs}.psort.bam'),
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 + 4096 * attempt,
    params:
        prefix = lambda wildcards, input: input.index.rsplit('.', 1)[0],
    shell:
        'bwa mem -t {threads} -R "@RG\\tID:{wildcards.library_id}\\tPL:Illumina\\tSM:{wildcards.sample}" '
            '{params.prefix} {input.reads} | '
            'samtools sort -m 1024M -@ 2 -O BAM > {output.tmp_bam} '
            ' && '
            'sambamba markdup -t {threads} --overflow-list-size 1000000 --io-buffer-size=512 {output.tmp_bam} {output.bam} '
            ' && '
            'samtools index {output.bam}'


rule collect_aligned_strandseq:
    input:
        lambda wildcards: expand(
            'output/alignments/sseq_to_{{graph_type}}_graph/{{sample}}/{library_id}_MAP-TO_{{graph}}.{{tigs}}.psort.mdup.bam{ext}',
            library_id=SSEQ_SAMPLE_TO_LIBS[wildcards.sample],
            ext=['', '.bai']
        )
    output:
        'output/alignments/sseq_to_{graph_type}_graph/{sample}_MAP-TO_{graph}.{tigs}.fofn',
    run:
        assert len(input) > 1
        with open(output[0], 'w') as dump:
            _ = dump.write('\n'.join(sorted(input)) + '\n')


rule write_haploclust_config_json:
    input:
        bam_aln = 'output/alignments/sseq_to_{graph_type}_graph/{sample}_MAP-TO_{graph}.{tigs}.fofn',
        sources = 'repos/haploclust/R'
    output:
        'output/haploclust/{sample}/sseq_to_{graph_type}_graph/{sample}_MAP-TO_{graph}.{tigs}/haploclust.cfg.json'
    params:
        threads = config['num_cpu_high'],
        input_type = 'bam',
        num_clusters = 80,
        em_iter = 2,
        num_alignments = 30000,
    run:
        import json

        hc_config = {
            'num_cpu': params.threads,
            'input_type': params.input_type,
            'num_clusters': params.num_clusters,
            'em_iter': params.em_iter,
            'num_alignments': params.num_alignments,
            'r_source_folder': input.sources  # special parameter
        }

        hc_config['output_folder'] = str(pathlib.Path(output[0]).parent / 'clustering')
        input_bam_files = []
        with open(input.bam_aln, 'r') as fofn:
            for line in fofn:
                if line.strip().endswith('.bai'):
                    continue
                input_bam_files.append(line.strip())
        hc_config['input_alignment_files'] = sorted(input_bam_files)
        hc_config['input_folder'] = str(pathlib.Path(input_bam_files[0]).parent)

        output_files = [
            'hard_clusters.RData',
            'soft_clusters.RData',
            'MLclust.data',
            'ss_clusters.data',
            'cluster_partners.txt',
            'wc_cells_clusters.data'
        ]
        for out in output_files:
            out_path = str(pathlib.Path(hc_config['output_folder']) / out)
            cfg_key = out.rsplit('.', 1)[0]
            hc_config[cfg_key] = out_path
        
        with open(output[0], 'w') as dump:
            _ = json.dump(hc_config, dump, ensure_ascii=True, indent=2, sort_keys=True)

    # END OF RUN BLOCK


SSEQ_CLUSTER=[f'cluster{i}' for i in range(1,24)]


rule run_haploclust_script:
    input:
        bam = lambda wildcards: expand(
            'output/alignments/sseq_to_{{graph_type}}_graph/{{sample}}/{library_id}_MAP-TO_{{graph}}.{{tigs}}.psort.mdup.bam',
            library_id=SSEQ_SAMPLE_TO_LIBS[wildcards.sample]
        ),
        bai = lambda wildcards: expand(
            'output/alignments/sseq_to_{{graph_type}}_graph/{{sample}}/{library_id}_MAP-TO_{{graph}}.{{tigs}}.psort.mdup.bam.bai',
            library_id=SSEQ_SAMPLE_TO_LIBS[wildcards.sample]
        ),
        cfg = 'output/haploclust/{sample}/sseq_to_{graph_type}_graph/{sample}_MAP-TO_{graph}.{tigs}/haploclust.cfg.json'
    output:
        ss_clust='output/haploclust/{sample}/sseq_to_{graph_type}_graph/{sample}_MAP-TO_{graph}.{tigs}/clustering/ss_clusters.data',
        hard_clust='output/haploclust/{sample}/sseq_to_{graph_type}_graph/{sample}_MAP-TO_{graph}.{tigs}/clustering/hard_clusters.RData',
        soft_clust='output/haploclust/{sample}/sseq_to_{graph_type}_graph/{sample}_MAP-TO_{graph}.{tigs}/clustering/soft_clusters.RData',
        ML_clust='output/haploclust/{sample}/sseq_to_{graph_type}_graph/{sample}_MAP-TO_{graph}.{tigs}/clustering/MLclust.data',
        ss_clust_sp=expand(
            'output/haploclust/{{sample}}/sseq_to_{{graph_type}}_graph/{{sample}}_MAP-TO_{{graph}}.{{tigs}}/clustering/ss_clusters_{cluster}.data',
            cluster=SSEQ_CLUSTER),
        clust_pairs='output/haploclust/{sample}/sseq_to_{graph_type}_graph/{sample}_MAP-TO_{graph}.{tigs}/clustering/clust_partners.txt',
        wc_cells_clusters='output/haploclust/{sample}/sseq_to_{graph_type}_graph/{sample}_MAP-TO_{graph}.{tigs}/clustering/wc_cells_clusters.data'
    log:
        'log/output/haploclust/{sample}/sseq_to_{graph_type}_graph/{sample}_MAP-TO_{graph}.{tigs}.log'
    benchmark:
        'rsrc/output/haploclust/{sample}/sseq_to_{graph_type}_graph/{sample}_MAP-TO_{graph}.{tigs}.rsrc'
    conda: '../../environment/conda/conda_rhaploclust.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: 12 * attempt,
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt
    params:
        input_type='bam',
        num_clusters=80,
        num_alignments=30000,
        EMiter=2,
    threads: config['num_cpu_high']
    script: str(pathlib.Path(pathlib.Path.cwd(), 'repos', 'haploclust', 'pipeline', 'utils', 'SaaRclust.snakemake.R'))
