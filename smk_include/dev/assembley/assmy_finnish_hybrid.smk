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
        assert 'MQ0' in wildcards.tigs
        if 'MQ0' in wildcards.tigs:
            mapq = 'mq00'
        else:
            raise
        if wildcards.tigs.endswith('XY'):
            chrom = 'chrXY'
        elif wildcards.tigs.endswith('X'):
            chrom = 'chrX'
        elif wildcards.tigs.endswith('Y'):
            chrom = 'chrY'
        else:
            raise ValueError(f'Unknown chrom: {str(wildcards)}')
        template = 'output/read_subsets/{chrom}/{sample_info}_{sample}_ONTUL.{chrom}-reads.{mapq}.fasta.gz'
        formatter = {
            'sample': wildcards.sample,
            'sample_info': wildcards.sample_info,
            'mapq': mapq,
            'chrom': chrom
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
    elif 'MQ0Y' in wildcards.tigs or 'MQ0XY' in wildcards.tigs:
        if any(x in wildcards.sample_info for x in ['DUO', 'TRIO', 'MIX']) or wildcards.sample.startswith('HC'):
            resources = config['num_cpu_high'], 262144, 23
        else:
            resources = config['num_cpu_high'], 65536, 11
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
        sample = CONSTRAINT_ALL_SAMPLES,
        ont_type = '(ONTUL|ONTEC|ONTHY)'
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
        runtime_hrs = lambda wildcards, attempt: attempt * attempt * attempt,
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt * attempt
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


def select_graph_to_dump(wildcards):
    assemblers = {
        'HAS': 'hifiasm',
        'MBG': 'mbg',
        'LJA': 'lja'
    }
    graph_paths = {
        ('hifiasm', 'input'): 'output/target_assembly/{chrom}/hifiasm/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.{tigs}.gfa',
        ('hifiasm', 'hybrid'): 'output/hybrid/110_final_graph/{sample_info}_{sample}.{ont_type}.{tigs}.final.gfa',
        ('mbg', 'input'): 'output/target_assembly/{chrom}/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.k{kmer}-w{window}-r{resolvek}.gfa',
        ('mbg', 'hybrid'): 'output/hybrid/110_final_graph/{sample_info}_{sample}.{ont_type}.{tigs}.final.gfa',
    }
    if wildcards.ont_type == 'UNSET':
        variant = 'input'
    elif wildcards.ont_type == 'ONTUL':
        variant = 'hybrid'
    else:
        raise ValueError(str(wildcards))

    assert 'MQ0Y' in wildcards.tigs or 'MQ0XY' in wildcards.tigs, f'Can only dump target assembly graphs: {wildcards.tigs}'
    mapq = 'mq00'

    formatter = dict(wildcards)
    formatter['mapq'] = mapq

    assembler = assemblers[wildcards.tigs.split('-')[0]]
    if variant == 'input':
        if 'OHEC' in wildcards.tigs:
            read_type = 'OHEC'
        elif 'OEC' in wildcards.tigs:
            read_type = 'ONTEC'
        elif 'HEC' in wildcards.tigs:
            read_type = 'HIFIEC'
        elif 'HAF' in wildcards.tigs:
            read_type = 'HIFIAF'
        else:
            raise
        formatter['read_type'] = read_type

        if wildcards.tigs.endswith('XY'):
            chrom = 'chrXY'
        elif wildcards.tigs.endswith('X'):
            chrom = 'chrX'
        elif wildcards.tigs.endswith('Y'):
            chrom = 'chrY'
        else:
            raise ValueError(f'Unknown chrom: {str(wildcards)}')
        formatter['chrom'] = chrom

        if assembler == 'hifiasm':
            formatter['tigs'] = get_hifiasm_tigs(wildcards.tigs)
        elif assembler == 'mbg':
            k, w, r = get_mbg_param(wildcards.tigs)
            formatter['kmer'] = k
            formatter['window'] = w
            formatter['resolvek'] = r
        else:
            raise
    
    graph_file_path = graph_paths[(assembler, variant)]
    gfa = graph_file_path.format(**formatter)
    return gfa


rule dump_final_graph_to_fasta:
    input:
        gfa = select_graph_to_dump
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


def select_fasta_to_align(wildcards):

    default_fasta = 'output/hybrid/200_final_post/{sample_info}_{sample}.{ont_type}.{tigs}.final.fasta'
    lja_fasta = 'output/target_assembly/chry_reads/lja/{sample_info}_{sample}_{read_type}.{mapq}.k{kmer}-K{resolvek}/assembly.fasta'
    formatter = dict(wildcards)
    if 'LJA' in wildcards.tigs:
        if 'MQ0' in wildcards.tigs:
            formatter['mapq'] = 'mq00'
        else:
            raise
        if 'OHEC' in wildcards.tigs:
            formatter['read_type'] = 'OHEC'
        elif 'OEC' in wildcards.tigs:
            formatter['read_type'] = 'ONTEC'
        elif 'HEC' in wildcards.tigs:
            formatter['read_type'] = 'HIFIEC'
        elif 'HAF' in wildcards.tigs:
            formatter['read_type'] = 'HIFIAF'
        else:
            raise
        small_k, large_k = get_lja_param(wildcards.tigs)
        formatter['kmer'] = small_k
        formatter['resolvek'] = large_k
        use_fasta = lja_fasta.format(**formatter)
    else:
        use_fasta = default_fasta.format(**formatter)
    return use_fasta


rule align_contigs_to_reference:
    input:
        fa_tigs = select_fasta_to_align,
        fa_ref = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta')
    output:
        paf = 'output/hybrid/210_align_ref/{sample_info}_{sample}_{ont_type}_{tigs}_MAP-TO_{reference}.paf.gz',
    benchmark:
        'rsrc/output/hybrid/210_align_ref/{sample_info}_{sample}_{ont_type}_{tigs}_MAP-TO_{reference}.mmap.rsrc',
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


def compute_overlap(row, ref_aln_start, ref_aln_end):    
    ovl = min(ref_aln_end, row['end']) - max(ref_aln_start, row['start'])
    return ovl


def compute_chry_seq_class(row, seq_classes):
    
    select_start = row['ref_aln_start'] < seq_classes['end']
    select_end = row['ref_aln_end'] > seq_classes['start']
    hit_classes = seq_classes.loc[select_start & select_end, :].copy()
    if hit_classes.empty:
        print(row)
        print(seq_classes)
        raise
    elif hit_classes.shape[0] > 1:
        hit_classes['overlap'] = hit_classes.apply(
            compute_overlap,
            axis=1,
            args=(row['ref_aln_start'], row['ref_aln_end'])
        )
        max_overlap = hit_classes['overlap'].max()
        hit_classes = hit_classes.loc[(hit_classes['overlap'] == max_overlap), :]
    assert hit_classes.shape[0] == 1
    df_idx = hit_classes.index[0]
    color = hit_classes.at[df_idx, 'color']
    seq_class = hit_classes.at[df_idx, 'class']
    order = hit_classes.at[df_idx, 'order']
    return color, seq_class, order


rule create_gfa_annotation:
    input:
        seq_classes = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.chrY-seq-classes.tsv'),
        aln_cache = 'output/hybrid/210_align_ref/{sample_info}_{sample}_{ont_type}_{tigs}_MAP-TO_{reference}.h5',
    output:
        csv = 'output/hybrid/220_gfa_annotation/{sample_info}_{sample}_{ont_type}_{tigs}_MAP-TO_{reference}.gfa-labels.csv',
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt
    run:
        import pandas as pd

        null_color = 'lightgrey'
        null_class = 'auto'
        null_order = 'any'

        seq_classes = pd.read_csv(
            input.seq_classes, sep='\t', header=None,
            names=['chrom', 'start', 'end', 'class', 'color', 'order'])
        seq_classes['end'] += 1

        keep_columns = [
            'qry_name',
            'qry_aln_fraction',
            'color',
            'seq_class',
            'seg_order',
            'ref_name'
        ]

        concat = []
        with pd.HDFStore(input.aln_cache, 'r') as hdf:
            for chrom in hdf.keys():
                ctg_aln = hdf[chrom]
                ctg_aln['qry_aln_length'] = ctg_aln['qry_aln_end'] - ctg_aln['qry_aln_start']
                ctg_aln['qry_aln_fraction'] = ctg_aln['qry_aln_length'] / ctg_aln['qry_length']
                if chrom.strip('/') != 'chrY':
                    ctg_aln['color'] = null_color
                    if chrom.strip('/') == 'chrX':
                        ctg_aln['seq_class'] = 'chrX'
                    else:
                        ctg_aln['seq_class'] = null_class
                    ctg_aln['seg_order'] = null_order
                else:
                    md = []
                    md_idx = []
                    for idx, row in ctg_aln.iterrows():
                        md.append(compute_chry_seq_class(row, seq_classes))
                        md_idx.append(idx)
                    md = pd.DataFrame.from_records(
                        md,
                        index=md_idx,
                        columns=['color', 'seq_class', 'seg_order']
                    )
                    ctg_aln = ctg_aln.merge(md, left_index=True, right_index=True)
                concat.append(ctg_aln[keep_columns])

        concat = pd.concat(concat, axis=0, ignore_index=False)
        concat.reset_index(drop=True, inplace=True)

        # merge split-align, all same region
        concat.drop_duplicates(
            ['qry_name', 'seg_order', 'ref_name'],
            keep='first',
            inplace=True
        )
        # drop split-align, keep max aligned
        concat.sort_values(['qry_name', 'qry_aln_fraction'], inplace=True, ascending=False)
        concat.drop_duplicates(
            ['qry_name'],
            keep='first',  # is sorted = keep max qry_aln_frac
            inplace=True
        )
        with open(output.csv, 'w') as dump:
            _ = dump.write(','.join(['Name', 'Chrom', 'AlnFrac', 'SeqClass', 'SegOrder', 'Color']) + '\n')
            concat[['qry_name', 'ref_name', 'qry_aln_fraction', 'seq_class', 'seg_order', 'color']].to_csv(
                dump,
                header=False,
                index=False
            )
        # END OF RUN BLOCK