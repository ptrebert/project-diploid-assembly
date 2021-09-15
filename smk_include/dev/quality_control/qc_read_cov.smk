
def set_alignment_memory(wildcards, attempt):
    """
    """
    base_mem = 24768
    if 'ONTEC' in wildcards.read_type:
        base_mem += 12288
    elif 'ONTUL' in wildcards.read_type:
        base_mem += 12288
    elif 'HIFIEC' in wildcards.read_type:
        pass
    elif 'HIFIAF' in wildcards.read_type:
        pass
    else:
        raise ValueError(str(wildcards))
    return base_mem * attempt


def set_alignment_runtime(wildcards, attempt):
    """
    NB: minimap2 w/ PAF output is very fast;
    unclear if the slowdown was always due to
    the piped SAM/BAM operations
    """
    base_hrs = 1
    if 'ONTEC' in wildcards.read_type:
        base_hrs += 2
    elif 'ONTUL' in wildcards.read_type:
        base_hrs += 2
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
        preset = 'map-hifi'
    elif 'ONTUL' in wildcards.read_type:
        preset = 'map-ont -k17'
    elif 'HIFIEC' in wildcards.read_type:
        preset = 'map-hifi'
    elif 'HIFIAF' in wildcards.read_type:
        preset = 'map-hifi'
    else:
        raise ValueError(str(wildcards))
    assert preset is not None
    return preset


rule qc_mmap_align_readsets:
    """
    https://github.com/lh3/minimap2/issues/771
    Above github issue contains some hints how to speed up alignment
    for ONT reads. Following this, set...
    -k17 [for ONTUL]
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
    shell:
        'minimap2 -t {resources.align_threads} -x {params.preset} --secondary=no '
        '--cap-kalloc=1g -K4g '
        '{input.reference} {input.reads} | pigz --best -p {resources.compress_threads} > {output}'


PAF_HEADER = [
    'read_name',
    'read_length',
    'read_aln_start',
    'read_aln_end',
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
        paf = 'output/alignments/reads_to_linear_ref/{sample}_{read_type}_{readset}_MAP-TO_{reference}.mmap.paf.gz',
        read_stats = 'input/{read_type}/{sample}_{read_type}_{readset}.stats.tsv.gz'
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

        metadata = dict()
        stats = pd.read_csv(
            input.read_stats,
            sep='\t',
            header=None,
            names=SEQTK_STATS_HEADER,
            usecols=[0,2,3,4,5,9]
        )
        metadata['total_reads'] = stats.shape[0]

        df = pd.read_csv(
            input.paf,
            sep='\t',
            header=None,
            names=PAF_HEADER,
            usecols=[0,1,2,3,4,5,7,8,9,10,11,16]
        )
        metadata['total_alignments'] = df.shape[0]
        metadata['total_aligned_reads'] = df['read_name'].nunique()

        unaligned = stats.loc[~stats['read_name'].isin(df['read_name']), :].copy()
        metadata['unaligned_reads'] = unaligned.shape[0]

        stats = stats.loc[~stats['read_name'].isin(unaligned['read_name']), :].copy()

        df = df.merge(stats, left_on='read_name', right_on='read_name', how='outer')
        df['read_aligned_length'] = (df['read_aln_end'] - df['read_aln_start']).astype('int64')
        df['read_aligned_pct'] = (df['read_aligned_length'] / df['read_length'] * 100).round(2)
        df['ref_aligned_length'] = (df['ref_aln_end'] - df['ref_aln_start']).astype('int64')
        df['divergence'] = df['divergence'].apply(lambda x: float(x.split(':')[-1]))
        assert df.notna().all(axis=0).all()

        with pd.HDFStore(output.cache, 'w', complevel=9, complib='blosc') as hdf:
            hdf.put(f'reference/{wildcards.reference}', pd.Series(chroms, dtype='int64'), format='fixed')
            hdf.put('metadata', pd.Series(metadata, dtype='int64'), format='fixed')
            if not unaligned.empty:
                hdf.put('unaligned', unaligned, format='fixed')
            for chrom, alignments in df.groupby('ref_name'):
                hdf.put(chrom, alignments, format='fixed')
    # END OF RUN BLOCK