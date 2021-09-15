rule meryl_query_only_kmer_db:
    """
    Create DB containing k-mers unique to the query sequences
    (the sequences for which the QV estimate should be computed)
    """
    input:
        query_db = 'output/kmer_smp_db/{sample}_{readset1}.meryl',
        reference_db = 'output/kmer_smp_db/{sample}_{readset2}.meryl'
    output:
        query_only = directory('output/kmer_op_db/{sample}_{readset1}_DIFF_{readset2}.meryl')
    benchmark:
        'rsrc/output/kmer_op_db/{sample}_{readset1}_DIFF_{readset2}.meryl.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * (attempt ** 4),
    shell:
        'meryl difference output {output.query_only} {input.query_db} {input.reference_db}'


def set_meryl_memory(wildcards, attempt):

    memory_mb = 8192 + 2048 * attempt
    if wildcards.read_type == 'ONTUL':
        # raw ONTUL reads
        memory_mb = 262144 + 65536 * attempt            
    return memory_mb


rule meryl_generate_individual_kmer_stats:
    """
    -existence:
        Generate a tab-delimited line for each input sequence with the
        number of kmers in the sequence, in the database and common to both.
    """
    input:
        query_only = 'output/kmer_op_db/{sample}_{read_type}_{readset}_DIFF_{readset2}.meryl',
        query_reads = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.read_type]
    output:
        table = 'output/kmer_stats/{sample}_{read_type}_{readset}_DIFF_{readset2}.seqkm.tsv'
    benchmark:
        'rsrc/output/kmer_stats/{sample}_{read_type}_{readset}_DIFF_{readset2}.seqkm.meryl.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: set_meryl_memory(wildcards, attempt),
        runtime_hrs = lambda wildcards, attempt: 8 * attempt
    shell:
        'meryl-lookup -existence -sequence {input.query_reads} -mers {input.query_only} > {output.table}'


def prob_base_correct(kmer_shared, kmer_total, kmer_size):
    return (kmer_shared / kmer_total) ** (1/kmer_size)


def base_error_rate(kmer_assembly_only, kmer_total, kmer_size):
    return 1 - (1 - kmer_assembly_only / kmer_total) ** (1/kmer_size)


def qv_estimate(error_rate):
    return -10 * math.log10(error_rate)


rule compute_global_query_qv_estimate:
    """
    Formulas for QV estimation as stated in 

    Rhie, A., Walenz, B.P., Koren, S. et al.
    Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies.
    Genome Biol 21, 245 (2020). https://doi.org/10.1186/s13059-020-02134-9

    Methods section "Consensus quality (QV) estimation"
    """
    input:
        query_stats = 'output/kmer_stats/{sample}_{readset1}.meryl.stats.h5',
        query_only_stats = 'output/kmer_stats/{sample}_{readset1}_DIFF_{readset2}.meryl.stats.h5'
    output:
        'output/qv_estimate/{sample}_{readset1}_REF_{readset2}.qv.tsv'
    run:
        import pandas as pd
        import math

        with pd.HDFStore(input.query_stats, mode='r') as hdf:
            num_kmer_query_total = hdf['statistics']['present']
            kmer_size_query_total = hdf['statistics']['kmer_size']

        with pd.HDFStore(input.query_only_stats, mode='r') as hdf:
            num_kmer_query_only = hdf['statistics']['present']
            kmer_size_query_only = hdf['statistics']['kmer_size']
        if kmer_size_query_total != kmer_size_query_only:
            raise ValueError(f'k-mer sizes do not match: {kmer_size_query_total} vs {kmer_size_query_only}')

        error_rate = base_error_rate(num_kmer_query_only, num_kmer_query_total, kmer_size_query_total)
        qv_est = round(-10 * math.log10(error_rate), 1)

        with open(output[0], 'w') as table:
            _ = table.write(f'sample\t{wildcards.sample}\n')
            _ = table.write(f'query_sequences\t{wildcards.readset1}\n')
            _ = table.write(f'reference_sequences\t{wildcards.readset2}\n')
            _ = table.write(f'kmer_size\t{kmer_size_query_only}\n')
            _ = table.write(f'kmer_query_only\t{num_kmer_query_only}\n')
            _ = table.write(f'kmer_query_total\t{num_kmer_query_total}\n')
            _ = table.write(f'error_rate\t{error_rate}\n')
            _ = table.write(f'QV_estimate\t{qv_est}\n')
    # END OF RUN BLOCK


rule compute_local_query_qv_estimate:
    """
    Comment regarding QV computation:
    The input meryl DB used for the lookup/existence operation contains
    only kmers unique to the (long read) sequences (i.e., the kmers from the
    short read dataset were removed [op difference]). Hence, the QV computation
    below uses "kmers shared between read sequence and read DB" as sequence-unique
    kmers (= not supported by Illumina/short read kmers), divided by the number
    of total kmers in the read sequence.
    """
    input:
        tsv = 'output/kmer_stats/{sample}_{readset}_DIFF_{ref_reads}.seqkm.tsv'
    output:
        hdf = 'output/qv_estimate/{sample}_{readset}_REF_{ref_reads}.seq-qv.h5'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    params:
        kmer_size = 31
    run:
        import pandas as pd
        import numpy as np

        table_header = ['read_name', 'kmer_in_seq', 'kmer_in_db', 'kmer_shared']
        df = pd.read_csv(input.tsv, sep='\t', names=table_header, header=None, index_col=None)

        # uncertain: for a test case of HiFi reads, the number of shared (= unsupported)
        # kmers was almost always 0. This might be an artifact due to the length
        # (statistical: expected number of errors per read for error rates ~10e-4 ~ 0),
        # or indicates some bug/problem...
        # For now, deal with that situation explicitly:

        df['error_rate'] = 0
        df['qv_estimate'] = -1

        select_nz = np.array(df['kmer_shared'] > 0, dtype=np.bool)

        df.loc[select_nz, 'error_rate'] = 1 - (1 - df.loc[select_nz, 'kmer_shared'] / df.loc[select_nz, 'kmer_in_seq']) ** (1/params.kmer_size)
        df.loc[select_nz, 'qv_estimate'] = (-10 * np.log10(df.loc[select_nz, 'error_rate'])).round(1)

        with pd.HDFStore(output.hdf, 'w', complevel=6) as hdf:
            hdf.put('sequence_qv', df, format='fixed')
