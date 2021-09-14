
def set_alignment_memory(wildcards, attempt):
    """
    aligning the ONT-UL reads is driving me crazy...
    now use the hammer
    """
    base_mem = 94208
    if 'ONTEC' in wildcards.read_type:
        base_mem += 24768
    elif 'ONTUL' in wildcards.read_type:
        base_mem = 262144
    elif 'HIFIEC' in wildcards.read_type:
        pass
    elif 'HIFIAF' in wildcards.read_type:
        pass
    else:
        raise ValueError(str(wildcards))
    return base_mem * attempt


def set_alignment_runtime(wildcards, attempt):
    """
    aligning the ONT-UL reads is driving me crazy...
    now use the hammer
    """
    base_hrs = 36
    if 'ONTEC' in wildcards.read_type:
        base_hrs = 72
    elif 'ONTUL' in wildcards.read_type:
        base_hrs = 167
    elif 'HIFIEC' in wildcards.read_type:
        pass
    elif 'HIFIAF' in wildcards.read_type:
        pass
    else:
        raise ValueError(str(wildcards))
    runtime_limit = base_hrs * attempt
    if runtime_limit > 167:
        raise ValueError(f'Runtime limit exceeds longq limit: {runtime_limit} / {str(wildcards)}')
    return runtime_limit


def set_alignment_preset(wildcards):
    preset = None
    if 'ONTEC' in wildcards.read_type:
        preset = 'map-pb'
    elif 'ONTUL' in wildcards.read_type:
        preset = 'map-ont'
    elif 'HIFIEC' in wildcards.read_type:
        preset = 'map-pb'
    elif 'HIFIAF' in wildcards.read_type:
        preset = 'map-pb'
    else:
        raise ValueError(str(wildcards))
    assert preset is not None
    return preset


rule qc_mmap_align_readsets:
    """
    https://github.com/lh3/minimap2/issues/771
    Above github issue contains some hints how to speed up alignment
    for ONT reads. Following this, set...
    -k17
    --cap-kalloc=1g
    -K4g
    """
    input:
        reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.read_type],
        reference = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta'),
    output:
        bam = 'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.mmap.paf.gz',
    log:
        'log/output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.mmap.log'
    benchmark:
        'rsrc/output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.mmap.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: set_alignment_memory(wildcards, attempt),
        runtime_hrs = lambda wildcards, attempt: set_alignment_runtime(wildcards, attempt),
        align_threads = config['num_cpu_high'] - config['num_cpu_low'],
        compress_threads = config['num_cpu_low'],
    params:
        preset = lambda wildcards: set_alignment_preset(wildcards),
        #temp_prefix = lambda wildcards: f'temp/mmap/{wildcards.reference}/{wildcards.sample}/{wildcards.read_type}/{wildcards.readset}/tmp_stsort_',
        #temp_dir = lambda wildcards: f'temp/mmap/{wildcards.reference}/{wildcards.sample}/{wildcards.read_type}/{wildcards.readset}/',
        validate = lambda wildcards, input:validate_readset(wildcards.readset, input.reads)
    shell:
        'minimap2 -t {resources.align_threads} -x {params.preset} --secondary=no '
        '-k17 --cap-kalloc=1g -K4g '
        '{input.reference} {input.reads} | pigz --best -p {resources.compress_threads} > {output}'


PAF_HEADER = [
    'query_name',
    'query_length',
    'query_start',
    'query_end',
    'orientation',
    'target_name',
    'target_length',
    'target_start',
    'target_end',
    'res_matches',
    'block_length',
    'mapq',
    'aln_type',
    'tag_cm',
    'tag_s1',
    'tag_s2',
    'divergence',
    'tag_rl'
]


BED_HEADER = [
    'chrom',
    'start',
    'end',
    'read_name',
    'mapq',
    'align_orient'
]

SEQTK_STATS_HEADER = [
    'read_name',
    'read_length',
    'num_A',
    'num_C',
    'num_G',
    'num_T',
    'ambig_2',
    'ambig_3',
    'ambig_4',
    'num_CpG',
    'tv',
    'ts',
    'num_CpG_ts'
]

rule cache_read_alignments:
    input:
        faidx = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta.fai'),
        bed = 'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.mmap.paf.gz',
        read_stats = 'input/{sample}_{read_type}_{readset}.stats.tsv.gz'
    output:
        cache = 'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.cov.cache.h5'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    run:
        import pandas as pd

        chroms = dict()
        with open(input.faidx, 'r') as faidx:
            for line in faidx:
                if not line.strip():
                    continue
                c, l = line.split()[:2]
                chroms[c] = int(l)
        chroms['genome'] = sum(chroms.values())
        chrom_sizes = pd.Series(chroms, dtype='int64')

        stats = pd.read_csv(
            input.read_stats,
            sep='\t',
            header=None,
            names=SEQTK_STATS_HEADER,
            usecols=[0,2,3,4,5,9]
        )

        df = pd.read_csv(input.bed, sep='\t', header=None, names=BED_HEADER)

        unaligned = stats.loc[~stats['read_name'].isin(df['read_name']), :].copy()
        stats = stats.loc[~stats['read_name'].isin(unaligned['read_name']), :].copy()

        df = df.merge(stats, left_on='read_name', right_on='read_name', how='outer')
        df['ref_aligned_length'] = df['end'] - df['start']
        df['ref_qry_fraction'] = (df['read_length'] / df['ref_aligned_length']).round(5)
        assert df.notna().all(axis=0).all()

        with pd.HDFStore(output.cache, 'w', complevel=9, complib='blosc') as hdf:
            hdf.put(f'reference/{wildcards.reference}', chrom_sizes, format='fixed')
            hdf.put('unaligned', unaligned, format='fixed')
            for chrom, alignments in df.groupby('chrom'):
                hdf.put(chrom, alignments, format='fixed')