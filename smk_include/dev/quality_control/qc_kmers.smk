
rule meryl_count_reference_kmers:
    input:
        sequence = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta')
    output:
        kmer_db = directory('/gpfs/project/projects/medbioinf/data/references/{reference}.k15.db'),
        rep_kmer = '/gpfs/project/projects/medbioinf/data/references/{reference}.k15.rep-grt09998.txt'
    benchmark:
        'rsrc/output/kmer_db/{reference}.meryl-ref-kmer.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 12288 + 4096 * attempt,
        mem_total_gb = lambda wildcards, attempt: int((12288 + 4096 * attempt)/1024) - 2,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        "meryl count k=15 threads={threads} memory={resources.mem_total_gb} output {output.kmer_db} {input.sequence} "
            " && "
            "meryl print greater-than distinct=0.9998 {output.kmer_db} > {output.rep_kmer}"


rule meryl_count_kmers:
    """
    Note: do not count k-mers in HPC-space b/c k-mer
    counts will only be used for QV estimation at the moment.
    Unsupported HP-errors should show up as such.
    """
    input:
        sequence = lambda wildcards: SAMPLE_INFOS[wildcards.sample][wildcards.read_type]
    output:
        kmer_db = directory('output/kmer_smp_db/{sample}_{read_type}_{readset}.meryl'),
    benchmark:
        'rsrc/output/kmer_smp_db/{sample}_{read_type}_{readset}.meryl.meryl.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: count_kmer_memory(wildcards, attempt),
        mem_total_gb = lambda wildcards, attempt: count_kmer_memory(wildcards, attempt, 'gb'),
        runtime_hrs = lambda wildcards, attempt: count_kmer_runtime(wildcards, attempt)
    params:
        kmer_size = 31,
    shell:
        "meryl count k={params.kmer_size} threads={threads} memory={resources.mem_total_gb} output {output.kmer_db} {input.sequence}"



def select_meryl_database(wildcards):

    if 'DIFF' in wildcards.readset:
        p = 'output/kmer_op_db/{readset}.meryl'
    else:
        p = 'output/kmer_smp_db/{readset}.meryl'
    return p


rule meryl_dump_db_statistics:
    """
    This operation is not documented in the cli help
    """
    input:
        kmer_db = select_meryl_database
    output:
        stats = 'output/kmer_stats/{readset}.meryl.statistics'
    benchmark:
        'rsrc/output/kmer_stats/{readset}.meryl.stats.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        "meryl statistics {input.kmer_db} > {output.stats}"


rule store_meryl_db_statistics:
    """
    NB: unique = singletons, k-mers with abundance 1

    FORMAT:

    Found 1 command tree.
    Number of 31-mers that are:
        unique              698268532  (exactly one instance of the kmer is in the input)
        distinct           2560577684  (non-redundant kmer sequences in the input)
        present           70996266663  (...)
        missing   4611686015866810220  (non-redundant kmer sequences not in the input)

                number of   cumulative   cumulative     presence
                distinct     fraction     fraction   in dataset
    frequency        kmers     distinct        total       (1e-6)
    --------- ------------ ------------ ------------ ------------
            1    698268532       0.2727       0.0098     0.000014
            2     20636896       0.2808       0.0104     0.000028
            3      6183975       0.2832       0.0107     0.000042
            4      3214166       0.2844       0.0109     0.000056
            5      2299957       0.2853       0.0110     0.000070
    """
    input:
        stats = 'output/kmer_stats/{readset}.meryl.statistics'
    output:
        hdf = 'output/kmer_stats/{readset}.meryl.stats.h5'
    benchmark:
        'rsrc/output/kmer_stats/{readset}.meryl.hdf.rsrc'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    run:
        import pandas as pd
        db_statistics = ['unique', 'distinct', 'present', 'missing']
        db_stats = {}
        kmer_freqs = []
        freq_line = False

        with open(input.stats, 'r') as txt:
            for ln, line in enumerate(txt, start=1):
                if not line.strip():
                    continue
                elif freq_line:
                    columns = line.strip().split()
                    assert len(columns) == 5, f'Malformed frequency line {ln}: {line.strip()}'
                    kmer_freqs.append(
                        (
                            int(columns[0]),
                            int(columns[1]),
                            float(columns[2]),
                            float(columns[3]),
                            float(columns[4])
                        )
                    )
                elif line.strip().startswith('Number of'):
                    parts = line.strip().split()
                    kmer_size = int(parts[2].split('-')[0])
                    db_stats['kmer_size'] = kmer_size
                else:
                    if any(line.strip().startswith(x) for x in db_statistics):
                        parts = line.strip().split()
                        statistic = parts[0]
                        assert statistic in db_statistics
                        try:
                            value = int(parts[1])
                        except ValueError:
                            # this happens b/c of "distinct" appearing
                            # twice (also in table header)
                            if parts[1] == 'fraction':
                                continue
                            raise
                        db_stats[statistic] = value
                    elif line.strip().startswith('-------'):
                        freq_line = True
                        continue
                    else:  # e.g., table header
                        continue
        assert len(db_stats) == 5, f'Missing DB stats: {db_stats}'
        db_stats = pd.Series(db_stats, index=None, name='db_statistics')
        kmer_freqs = pd.DataFrame.from_records(
            kmer_freqs,
            columns=[
                'frequency',
                'num_distinct_kmers',
                'cum_fraction_distinct',
                'cum_fraction_total',
                'presence_in_dataset'
            ]
        )
        with pd.HDFStore(output.hdf, 'w', complevel=9) as hdf:
            hdf.put('statistics', db_stats, format='fixed')
            hdf.put('kmer_freq', kmer_freqs, format='fixed')
    # END OF RUN BLOCK