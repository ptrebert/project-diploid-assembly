import pathlib as pl

"""
Follow HiFi/ONT hybrid approach developed by MR
https://github.com/maickrau/hybrid-assembly/blob/master/commands.sh
[lines 46-65]
Around commit e932291, commands.sh was shortened the relevant steps are now in lines
14-35
"""


def select_hybrid_assm_ont_reads(wildcards):

    if 'MQ' in wildcards.tigs:
        assert wildcards.ont_type == 'ONTUL'
        assert 'YRAW' in wildcards.tigs
        if 'MQ0' in wildcards.tigs:
            mapq = 'mq00'
        else:
            raise
        template = 'output/read_subsets/chry/{sample_info}_{sample}_ONTUL.chrY-reads.{mapq}.fasta.gz'
        formatter = {
            'sample': wildcards.sample,
            'sample_info': wildcards.sample_info,
            'mapq': mapq
        }
        ont_reads = template.format(**formatter)
    else:
        ont_reads = str(SAMPLE_INFOS[wildcards.sample][wildcards.ont_type])
    if 'gpfs' in ont_reads:
        ont_reads = str(pl.Path('/hilbert', ont_reads.strip('/')))
    return ont_reads


def set_graphaligner_hybrid_resources(wildcards):

    if wildcards.tigs in ['TIGRAW', 'TIGPRI', 'TIGALT']:
        resources = config['num_cpu_max'], 303104, 167
    elif 'YRAW' in wildcards.tigs:
        resources = config['num_cpu_high'], 24576, 2
    else:
        raise
    return resources


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
        reads = select_hybrid_assm_ont_reads
    output:
        gaf = 'output/hybrid/ont_to_graph/{sample_info}_{sample}.{ont_type}.{tigs}.gaf',
        hybrid_reads = 'output/hybrid/ont_to_graph/{sample_info}_{sample}.{ont_type}.{tigs}.ONTHY.fasta',
    log:
        'log/output/hybrid/ont_to_graph/{sample_info}_{sample}.{ont_type}.{tigs}.ga.log',
    benchmark:
        'rsrc/output/hybrid/ont_to_graph/{sample_info}_{sample}.{ont_type}.{tigs}.ga.rsrc',
    wildcard_constraints:
        sample = CONSTRAINT_ALL_SAMPLES
#    conda: '../../../environment/conda/conda_biotools.yml'
    threads: lambda wildcards: set_graphaligner_hybrid_resources(wildcards)[0]
    resources:
        mem_total_mb = lambda wildcards, attempt: set_graphaligner_hybrid_resources(wildcards)[1] * attempt,
        runtime_hrs = lambda wildcards, attempt: set_graphaligner_hybrid_resources(wildcards)[2] * attempt
    shell:
        'module load Singularity && singularity exec '
        '--bind /:/hilbert {input.container} '
        'GraphAligner -g {input.graph} -f {input.reads} '
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
        raw_unitigs = 'output/hybrid/afr_mix/hifiasm/{afr_sample1}-U-{afr_sample2}.r_utg.gfa'
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


rule dump_afr_mix_graph_to_fasta:
    input:
        gfa = 'output/hybrid/afr_mix/hifiasm/{afr_sample1}-U-{afr_sample2}.r_utg.gfa'
    output:
        stats = 'output/hybrid/afr_mix/{afr_sample1}-U-{afr_sample2}.TIGRAW.gfa.stats.txt',
        fasta = 'output/hybrid/afr_mix/{afr_sample1}-U-{afr_sample2}.TIGRAW.fasta',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
    shell:
        'gfatools stat {input.gfa} > {output.stats}'
            ' && '
        'gfatools gfa2fa {input.gfa} > {output.fasta}'


rule align_afr_mix_contigs_to_reference:
    input:
        fa_tigs = 'output/hybrid/afr_mix/{afr_sample1}-U-{afr_sample2}.TIGRAW.fasta',
        fa_ref = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta')
    output:
        paf = 'output/hybrid/afr_mix/ctg_ref_align/{afr_sample1}-U-{afr_sample2}.TIGRAW_MAP-TO_{reference}.paf.gz',
    benchmark:
        'rsrc/output/hybrid/afr_mix/ctg_ref_align/{afr_sample1}-U-{afr_sample2}.TIGRAW_MAP-TO_{reference}.mmap.rsrc',
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 24576 + 16384 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
    shell:
        'minimap2 -t {threads} '
            '--cap-kalloc=1g -K4g '
            '--secondary=no -x asm20 -m 10000 --end-bonus=100 '
            '{input.fa_ref} {input.fa_tigs} | '
        'pigz -p 4 --best > {output.paf}'


rule cache_afr_mix_contig_to_reference_alignment:
    input:
        paf = 'output/hybrid/afr_mix/ctg_ref_align/{afr_sample1}-U-{afr_sample2}.TIGRAW_MAP-TO_{reference}.paf.gz',
    output:
        hdf = 'output/hybrid/afr_mix/ctg_ref_align/{afr_sample1}-U-{afr_sample2}.TIGRAW_MAP-TO_{reference}.h5',
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd
        import collections as col
        df = pd.read_csv(input.paf, sep='\t', header=None, names=PAF_HEADER, usecols=PAF_USE)
        df['divergence'] = (df['divergence'].apply(lambda x: float(x.split(':')[-1]))).astype('float16')
        orientation_map = col.defaultdict(int)
        orientation_map['+'] = 1
        orientation_map['-'] = -1
        df['orientation'] = df['orientation'].replace(orientation_map)
        df['orientation'] = df['orientation'].astype('int8')

        with pd.HDFStore(output.hdf, 'w', complevel=9, complib='blosc') as hdf:
            for ref_chrom, alignments in df.groupby('ref_name'):
                hdf.put(f'{ref_chrom}', alignments, format='fixed')
    # END OF RUN BLOCK


#######################################
### BELOW: hybrid assembly approach ###
#######################################


HYBRID_SCRIPT_PATH = 'repos/hybrid-assembly/scripts'


rule filter_ont_to_graph_alignment:
    input:
        gaf = 'output/hybrid/ont_to_graph/{sample_info}_{sample}.{ont_type}.{tigs}.gaf',
    output:
        gaf = 'output/hybrid/10_mapq_length_filter/{sample_info}_{sample}.{ont_type}.{tigs}.qlfilter.gaf',
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt,
    params:
        min_aligned_length = 0.8,
        min_mapq = 20
    shell:
        "awk -F '\\t' '{{if ($4-$3 >= $2*{params.min_aligned_length} && $12 >= {params.min_mapq}) print;}}' < {input.gaf} > {output.gaf}"


rule trim_graph_alignment:
    input:
        gaf = 'output/hybrid/10_mapq_length_filter/{sample_info}_{sample}.{ont_type}.{tigs}.qlfilter.gaf',
        graph = 'output/clean_graphs/{sample_info}_{sample}.{tigs}.gfa',
    output:
        gaf = 'output/hybrid/20_trim_graph_alignment/{sample_info}_{sample}.{ont_type}.{tigs}.trimmed.gaf'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = os.path.join(HYBRID_SCRIPT_PATH, 'trim_dbg_alignment.py'),
        edge_trim = 1500
    shell:
        '{params.script_exec} {input.graph} {params.edge_trim} < {input.gaf} > {output.gaf}'


rule calculate_node_coverage:
    input:
        graph = 'output/clean_graphs/{sample_info}_{sample}.{tigs}.gfa',
        trimmed_aln = 'output/hybrid/20_trim_graph_alignment/{sample_info}_{sample}.{ont_type}.{tigs}.trimmed.gaf'
    output:
        table = 'output/hybrid/30_node_coverages/{sample_info}_{sample}.{ont_type}.{tigs}.nodecov.csv'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = os.path.join(HYBRID_SCRIPT_PATH, 'calculate_coverage.py'),
        edge_trim = 1500
    shell:
        '{params.script_exec} {input.graph} < {input.trimmed_aln} > {output.table}'


rule extract_qlfilter_trimmed_paths:
    input:
        gaf = 'output/hybrid/20_trim_graph_alignment/{sample_info}_{sample}.{ont_type}.{tigs}.trimmed.gaf'
    output:
        listing = 'output/hybrid/40_qlfilter_trimmed_paths/{sample_info}_{sample}.{ont_type}.{tigs}.paths.txt'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'cut -f 6 < {input.gaf} > {output.listing}'


rule mapq_filter_ont_graph_alignment:
    input:
        gaf = 'output/hybrid/ont_to_graph/{sample_info}_{sample}.{ont_type}.{tigs}.gaf',
    output:
        gaf = 'output/hybrid/10_mapq_only_filter/{sample_info}_{sample}.{ont_type}.{tigs}.mqfilter.gaf',
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt,
    params:
        min_mapq = 20
    shell:
        "awk -F '\\t' '{{if ($12 >= {params.min_mapq}) print;}}' < {input.gaf} > {output.gaf}"


rule insert_alignment_gaps:
    input:
        gaf = 'output/hybrid/10_mapq_only_filter/{sample_info}_{sample}.{ont_type}.{tigs}.mqfilter.gaf',
        graph = 'output/clean_graphs/{sample_info}_{sample}.{tigs}.gfa',
    output:
        graph = 'output/hybrid/50_insert_aln_gaps/{sample_info}_{sample}.{ont_type}.{tigs}.gapped.gfa'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    params:
        script_exec = os.path.join(HYBRID_SCRIPT_PATH, 'insert_aln_gaps.py'),
        min_gap_coverage = 3,
        max_end_clip = 50
    shell:
        '{params.script_exec} {input.graph} {params.min_gap_coverage} {params.max_end_clip} < {input.gaf} > {output.graph}'


rule identify_unique_nodes:
    input:
        gapped_graph = 'output/hybrid/50_insert_aln_gaps/{sample_info}_{sample}.{ont_type}.{tigs}.gapped.gfa',
        trimmed_aln = 'output/hybrid/20_trim_graph_alignment/{sample_info}_{sample}.{ont_type}.{tigs}.trimmed.gaf'
    output:
        listing = 'output/hybrid/60_id_unique_nodes/{sample_info}_{sample}.{ont_type}.{tigs}.unique-nodes.txt'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = os.path.join(HYBRID_SCRIPT_PATH, 'estimate_unique_local.py'),
        long_node_threshold = 100000,
        solid_edge_threshold = 30,
        path_consistency_threshold = 0.8
    shell:
        '{params.script_exec} {input.gapped_graph} {input.trimmed_aln} '
            '{params.long_node_threshold} {params.solid_edge_threshold} {params.path_consistency_threshold} '
            '> {output.listing}'


rule find_bridges:
    input:
        uniq_nodes = 'output/hybrid/60_id_unique_nodes/{sample_info}_{sample}.{ont_type}.{tigs}.unique-nodes.txt',
        paths = 'output/hybrid/40_qlfilter_trimmed_paths/{sample_info}_{sample}.{ont_type}.{tigs}.paths.txt'
    output:
        listing = 'output/hybrid/70_bridges/{sample_info}_{sample}.{ont_type}.{tigs}.bridges.txt'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = os.path.join(HYBRID_SCRIPT_PATH, 'find_bridges.py'),
    shell:
        '{params.script_exec} {input.uniq_nodes} < {input.paths} > {output.listing}'


rule discard_invalid_bridges:
    input:
        bridges = 'output/hybrid/70_bridges/{sample_info}_{sample}.{ont_type}.{tigs}.bridges.txt'
    output:
        listing = 'output/hybrid/80_valid_bridges/{sample_info}_{sample}.{ont_type}.{tigs}.valid-bridges.txt'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = os.path.join(HYBRID_SCRIPT_PATH, 'remove_wrong_connections_2.py'),
    shell:
        'grep -v "(" {input.bridges} | grep -vP "^$" | {params.script_exec} | sort > {output.listing}'


rule select_majority_bridge:
    input:
        valid_bridges = 'output/hybrid/80_valid_bridges/{sample_info}_{sample}.{ont_type}.{tigs}.valid-bridges.txt'
    output:
        listing = 'output/hybrid/90_majority_bridges/{sample_info}_{sample}.{ont_type}.{tigs}.majority-bridges.txt'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = os.path.join(HYBRID_SCRIPT_PATH, 'pick_majority_bridge.py'),
    shell:
        '{params.script_exec} < {input.valid_bridges} > {output.listing}'


rule identify_forbidden_end:
    input:
        uniq_nodes = 'output/hybrid/60_id_unique_nodes/{sample_info}_{sample}.{ont_type}.{tigs}.unique-nodes.txt',
        gapped_graph = 'output/hybrid/50_insert_aln_gaps/{sample_info}_{sample}.{ont_type}.{tigs}.gapped.gfa',
        majority_bridges = 'output/hybrid/90_majority_bridges/{sample_info}_{sample}.{ont_type}.{tigs}.majority-bridges.txt',
        paths = 'output/hybrid/40_qlfilter_trimmed_paths/{sample_info}_{sample}.{ont_type}.{tigs}.paths.txt',
        nodecov = 'output/hybrid/30_node_coverages/{sample_info}_{sample}.{ont_type}.{tigs}.nodecov.csv',
    output:
        listing = 'output/hybrid/100_forbidden_end/{sample_info}_{sample}.{ont_type}.{tigs}.forbidden-ends.txt'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = os.path.join(HYBRID_SCRIPT_PATH, 'forbid_unbridged_tangles.py'),
        min_solid_coverage = 30
    shell:
        '{params.script_exec} {input.uniq_nodes} {input.gapped_graph} {input.majority_bridges} '
            '{input.paths} {input.nodecov} {params.min_solid_coverage} > {output.listing}'


rule build_connected_graph:
    input:
        gapped_graph = 'output/hybrid/50_insert_aln_gaps/{sample_info}_{sample}.{ont_type}.{tigs}.gapped.gfa',
        forbidden_ends = 'output/hybrid/100_forbidden_end/{sample_info}_{sample}.{ont_type}.{tigs}.forbidden-ends.txt',
        majority_bridges = 'output/hybrid/90_majority_bridges/{sample_info}_{sample}.{ont_type}.{tigs}.majority-bridges.txt',
    output:
        gfa = 'output/hybrid/110_final_graph/{sample_info}_{sample}.{ont_type}.{tigs}.final.gfa'
    log:
        'log/output/hybrid/110_final_graph/{sample_info}_{sample}.{ont_type}.{tigs}.build-final.log'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt
    params:
        script_exec = os.path.join(HYBRID_SCRIPT_PATH, 'connect_uniques.py'),
    shell:
        '{params.script_exec} {input.gapped_graph} {input.forbidden_ends} {input.majority_bridges} > {output.gfa} 2> {log}'


rule dump_final_graph_to_fasta:
    input:
        gfa = 'output/hybrid/110_final_graph/{sample_info}_{sample}.{ont_type}.{tigs}.final.gfa'
    output:
        stats = 'output/hybrid/200_final_post/{sample_info}_{sample}.{ont_type}.{tigs}.gfa.stats.txt',
        fasta = 'output/hybrid/200_final_post/{sample_info}_{sample}.{ont_type}.{tigs}.final.fasta',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
    shell:
        'gfatools stat {input.gfa} > {output.stats}'
            ' && '
        'gfatools gfa2fa {input.gfa} > {output.fasta}'


rule align_contigs_to_reference:
    input:
        fa_tigs = 'output/hybrid/200_final_post/{sample_info}_{sample}.{ont_type}.{tigs}.final.fasta',
        fa_ref = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta')
    output:
        paf = 'output/hybrid/210_align_ref/{sample_info}_{sample}_{ont_type}_{tigs}_MAP-TO_{reference}.paf.gz',
    benchmark:
        'rsrc/output/hybrid/210_align_ref/{sample_info}_{sample}_{ont_type}_{tigs}_MAP-TO_{reference}.mmap.rsrc',
    wildcard_constraints:
        tigs = '(TIGPRI|TIGALT|TIGRAW|ECMQ0YRAW|AFMQ0YRAW)'
    conda: '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 24576 + 16384 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
    shell:
        'minimap2 -t {threads} '
            '--cap-kalloc=1g -K4g '
            '--secondary=no -x asm20 -m 10000 --end-bonus=100 '
            '{input.fa_ref} {input.fa_tigs} | '
        'pigz -p 4 --best > {output.paf}'


PAF_HEADER = [
    'qry_name',
    'qry_length',
    'qry_aln_start',
    'qry_aln_end',
    'orientation',
    'ref_name',
    'ref_length',
    'ref_aln_start',
    'ref_aln_end',
    'residue_matches',
    'block_length',
    'mapq',
    'aln_type',
    'tag_cm',
    'tag_s1',
    'tag_s2',
    'divergence',
    'tag_rl'
]


PAF_USE = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16]


rule cache_contig_to_reference_alignment:
    input:
        paf = 'output/hybrid/210_align_ref/{sample_info}_{sample}_{ont_type}_{tigs}_MAP-TO_{reference}.paf.gz',
    output:
        hdf = 'output/hybrid/210_align_ref/{sample_info}_{sample}_{ont_type}_{tigs}_MAP-TO_{reference}.h5',
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt
    run:
        import pandas as pd
        import collections as col
        df = pd.read_csv(input.paf, sep='\t', header=None, names=PAF_HEADER, usecols=PAF_USE)
        df['divergence'] = (df['divergence'].apply(lambda x: float(x.split(':')[-1]))).astype('float16')
        orientation_map = col.defaultdict(int)
        orientation_map['+'] = 1
        orientation_map['-'] = -1
        df['orientation'] = df['orientation'].replace(orientation_map)
        df['orientation'] = df['orientation'].astype('int8')

        with pd.HDFStore(output.hdf, 'w', complevel=9, complib='blosc') as hdf:
            for ref_chrom, alignments in df.groupby('ref_name'):
                hdf.put(f'{ref_chrom}', alignments, format='fixed')
    # END OF RUN BLOCK
