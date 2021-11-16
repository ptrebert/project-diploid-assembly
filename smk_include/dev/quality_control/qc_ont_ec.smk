import pathlib as pl

def set_errormasking(wildcards):

    if 'HIFIAF' in wildcards.read_type:
        errm = 'collapse-msat'
    elif 'HIFIEC' in wildcards.read_type:
        errm = 'no'
    else:
        raise ValueError(str(wildcards))
    return errm


def get_read_path(sample, read_type, prefix=''):

    file_path = str(SAMPLE_INFOS[sample][read_type])
    if prefix:
        file_path = str(pl.Path(prefix, file_path.strip('/')))
        assert 'gpfs' in file_path, f'Not an abs path: {file_path}'
    return file_path


def set_mbg_resources(wildcards):
    """
    Update to MBG v1.0.7+, use default for all read types,
    unclear how much error-masking is affecting runtime
    for HIFI-AF read sets
    """
    if wildcards.read_type == 'HIFIAF':
        # uses error-masking, much faster apparently
        cpu, memory, runtime = config['num_cpu_high'], 114688, 23
    elif wildcards.read_type == 'HIFIEC':
        cpu, memory, runtime = config['num_cpu_high'], 385024, 72
    else:
        raise ValueError(str(wildcards))
    return cpu, memory, runtime


rule build_hifi_read_dbg:
    """
    sif = ancient('mbg.master.sif'),
    sif = ancient('mbg.UnitigResolve.sif'),

    'module load Singularity && singularity exec '
    '--bind /:/hilbert {input.sif} '

    NB: only works with MBG 1.0.7+

    """
    input:
        reads = lambda wildcards: get_read_path(wildcards.sample, wildcards.read_type)
    output:
        graph = 'output/mbg_hifi/{sample}_{read_type}_{readset}.MBG-k{kmer}-w{window}-r{resolve}.gfa',
        paths = 'output/mbg_hifi/{sample}_{read_type}_{readset}.MBG-k{kmer}-w{window}-r{resolve}.gaf'
    log:
        'log/output/mbg_hifi/{sample}_{read_type}_{readset}.MBG-k{kmer}-w{window}-r{resolve}.MBG.log'
    benchmark:
        'rsrc/output/mbg_hifi/{sample}_{read_type}_{readset}.MBG-k{kmer}-w{window}-r{resolve}.MBG.rsrc'
    wildcard_constraints:
        read_type = '(HIFIEC|HIFIAF)'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: lambda wildcards: set_mbg_resources(wildcards)[0]
    resources:
        mem_total_mb = lambda wildcards, attempt: set_mbg_resources(wildcards)[1] * attempt,
        runtime_hrs = lambda wildcards, attempt: min(set_mbg_resources(wildcards)[2] * attempt, 167)
    params:
        masking = set_errormasking,
        #input_path = lambda wildcards: get_read_path(wildcards.sample, wildcards.read_type, '/hilbert')
    shell:
        'MBG -i {input.reads} -t {threads} '
            '-k {wildcards.kmer} -w {wildcards.window} --resolve-maxk {wildcards.resolve} '
            '--error-masking {params.masking} '
            '--out {output.graph} --output-sequence-paths {output.paths} &> {log}'


def set_read_hpc(wildcards):

    if 'HIFIAF' in wildcards.graph_reads:
        hpc = '--hpc-collapse-reads'
    elif 'HIFIEC' in wildcards.graph_reads:
        hpc = ''
    else:
        raise ValueError(str(wildcards))
    return hpc


def set_graphaligner_resources(wildcards):

    # if int(wildcards.kmer) < 999:
    #     cpu, memory, runtime = config['num_cpu_max'], 786432, 167
    # else:
    cpu, memory, runtime = config['num_cpu_high'] + config['num_cpu_medium'], 114688, 23
    return cpu, memory, runtime


rule ont_to_graph_alignment:
    """
    """
    input:
        sif = ancient('graphaligner.MultiseedClusters.sif'),
        graph = 'output/mbg_hifi/{sample}_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}-r{resolve}.gfa',
        reads = 'input/ONTUL/{sample}_ONTUL_{readset}.fasta.gz'
    output:
        gaf = 'output/alignments/ont_to_mbg_graph/{sample}_ONTUL_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}-r{resolve}.gaf',
        ec_reads_clip = 'output/alignments/ont_to_mbg_graph/{sample}_ONTEC_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}-r{resolve}.fasta.gz',
        ec_reads_full = 'output/alignments/ont_to_mbg_graph/{sample}_ONTHY_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}-r{resolve}.fasta.gz'
    log:
        'log/output/alignments/ont_to_mbg_graph/{sample}_ONTUL_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}-r{resolve}.GA.log'
    benchmark:
        'rsrc/output/alignments/ont_to_mbg_graph/{sample}_ONTUL_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}-r{resolve}.GA.rsrc'
    threads: lambda wildcards: set_graphaligner_resources(wildcards)[0]
    resources:
        mem_total_mb = lambda wildcards, attempt: set_graphaligner_resources(wildcards)[1] * attempt,
        runtime_hrs = lambda wildcards, attempt: set_graphaligner_resources(wildcards)[2] * attempt
    params:
        preset = 'dbg',
        hpc = set_read_hpc,
    shell:
        'module load Singularity && singularity exec '
        '--bind /:/hilbert {input.sif} '
        'GraphAligner -g {input.graph} -f {input.reads} '
            '-x {params.preset} -t {threads} '
            '--min-alignment-score 5000 --multimap-score-fraction 1 '
            ' {params.hpc} '
            '--corrected-clipped-out {output.ec_reads_clip} '
            '--corrected-out {output.ec_reads_full} '
            '-a {output.gaf} &> {log}'


rule compute_ont_corrected_stats:
    input:
        'output/alignments/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}-r{resolve}.fasta.gz',
    output:
        'output/alignments/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}-r{resolve}.stats.tsv.gz',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: 4
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: 12 * attempt
    shell:
        'pigz -p 2 -d -c {input} | seqtk comp | pigz -p 2 --best > {output} '


rule cache_ont_corrected_read_stats:
    input:
        tsv = 'output/alignments/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}-r{resolve}.stats.tsv.gz'
    output:
        h5 = 'output/alignments/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}-r{resolve}.stats.h5'
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt
    run:
        import pandas as pd
        import pathlib as pl
        seqtk_columns = [
            'read_name',
            'read_length',
            'num_A',
            'num_C',
            'num_G',
            'num_T',
            'num_IUPAC2',
            'num_IUPAC3',
            'num_N',
            'num_CpG',
            'num_ts',
            'num_tv',
            'num_CpG_ts',
        ]
        input_file_name = pl.Path(input.tsv).name
        df = pd.read_csv(input.tsv, sep='\t', header=None, names=seqtk_columns)
        df['MBG_kmer'] = int(wildcards.kmer)
        df['MBG_window'] = int(wildcards.window)
        df['MBG_resolveK'] - int(wildcards.resolve)
        df['file'] = input_file_name
        with pd.HDFStore(output.h5, 'w', complib='blosc', complevel=9) as hdf:
            hdf.put('cache', df, format='fixed')
    # END OF RUN BLOCK


# rule deduplicate_ont_corrected_reads:
#     input:
#         cache = expand(
#             'output/alignments/ont_to_mbg_graph/{{sample}}_{{read_type}}_{{readset}}_MAP-TO_{{graph_reads}}_{{graph_readset}}.MBG-k{kmer}-w{window}.stats.h5',
#             zip,
#             kmer=config['mbg_kmers'],
#             window=config['mbg_windows']
#         )
#     output:
#         expand(
#             'output/alignments/ont_to_mbg_graph/{{sample}}_{{read_type}}_{{readset}}_MAP-TO_{{graph_reads}}_{{graph_readset}}.MBG-k{kmer}-w{window}.select-reads.txt',
#             zip,
#             kmer=config['mbg_kmers'],
#             window=config['mbg_windows']
#         )
#     resources:
#         mem_total_mb = lambda wildcards, attempt: 8192 * attempt
#     run:
#         import pandas as pd
#         import pathlib as pl
#         import collections as col

#         df = []
#         out_path = ''
#         for cache_file in input.cache:
#             out_path = pl.Path(cache_file).parent
#             df.append(pd.read_hdf(cache_file, 'cache'))
#         df = pd.concat(df, axis=0, ignore_index=False)
#         multiplicity = df['read_name'].value_counts()
#         df['multiplicty'] = df['read_name'].apply(lambda x: multiplicity[x])

#         out_files = dict()
#         for file_name, read_names in df.loc[df['multiplicity'] < 2, :].groupby('file')['read_name']:
#             try:
#                 out_file_name = out_files[file_name]
#             except KeyError:
#                 out_file_name = file_name.replace('.stats.tsv.gz', '.select-reads.txt')
#                 out_files[file_name] = out_file_name
#             with open(pl.Path(out_path, out_file_name), 'w') as dump:
#                 _ = dump.write('\n'.join(read_names.tolist()) + '\n')
        
#         out_cache = col.defaultdict(list)
#         for read_name, file_names in df.loc[df['multiplicity'] > 1, :].groupby('read_name')['file']:
#             select_file = file_names.sample(1).values[0]
#             out_cache[out_files[select_file]].append(read_name)

#         for out_file_name, read_names in out_cache.items():
#             with open(pl.Path(out_path, out_file_name), 'a') as dump:
#                 _ = dump.write('\n'.join(read_names.tolist()) + '\n')
        
#         # END OF RUN BLOCK


# rule extract_selected_ontec_reads:
#     input:
#         reads = 'output/alignments/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.fasta.gz',
#         names = 'output/alignments/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.select-reads.txt',
#     output:
#         'output/dedup_ontec/ont_to_mbg_graph/{sample}_{read_type}_{readset}_MAP-TO_{graph_reads}_{graph_readset}.MBG-k{kmer}-w{window}.uniq.fasta.gz',
#     conda:
#         '../../../environment/conda/conda_biotools.yml'
#     wildcard_constraints:
#         graph_reads = '(HIFIEC|HIFIAF)'
#     threads: config['num_cpu_low']
#     resources:
#         mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
#         runtime_hrs = lambda wildcards, attempt: 4 * attempt
#     shell:
#         'pigz -p 2 -d -c {input.reads} | seqtk subseq /dev/stdin {input.names} | pigz -p {threads} --best > {output}'


# rule merge_extracted_ontec_reads:
#     input:
#         subsets = expand(
#             'output/dedup_ontec/ont_to_mbg_graph/{{sample}}_{{read_type}}_{readset}_MAP-TO_{{graph_reads}}_{{graph_readset}}.MBG-k{kmer}-w{window}.uniq.fasta.gz',
#             zip,
#             kmer=config['mbg_kmers'],
#             window=config['mbg_windows'],
#             readset=['guppy-5.0.11-sup-prom'] * len(config['mbg_kmers']),
#         )
#     output:
#         'input/{read_type}/{sample}_{read_type}_{graph_reads}-{graph_readset}.fasta.gz'
#     conda:
#         '../../../environment/conda/conda_biotools.yml'
#     wildcard_constraints:
#         graph_reads = '(HIFIEC|HIFIAF)'
#     threads: config['num_cpu_low']
#     resources:
#         mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
#         runtime_hrs = lambda wildcards, attempt: 4 * attempt
#     shell:
#         'pigz -p 2 -d -c {input.subsets} | pigz -p {threads} --best > {output}'
