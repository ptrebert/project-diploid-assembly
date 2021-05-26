
REFERENCE_FOLDER = '/beeond/data/references'

MIN_CONTIG_SIZE = 20000
MIN_GAP_SIZE = 1000
MAX_GAP_SIZE = 200000


rule master:
    input:
        'output/assemblies/T2Tv11_randsplit.fasta.fai',
        'output/assemblies/T2Tv11_sdsplit.fasta.fai',
        'output/kmer_db/T2Tv11_randsplit.k15.no-hpc.rep-grt09998.txt.gz',
        'output/kmer_db/T2Tv11_sdsplit.k15.no-hpc.rep-grt09998.txt.gz'


def rand_split_sequence(seq_name, sequence, min_gap, max_gap, min_contig):

    import io
    import random as rand

    sequence_length = len(sequence)

    max_contig = sequence_length // 3
    split_buffer = io.StringIO()

    last_end = 0
    order_number = 0
    while 1:
        start = last_end + rand.randint(min_gap, max_gap)
        if start > sequence_length:
            break
        contig_size = rand.randint(min_contig, max_contig)
        if start + contig_size > sequence_length:
            end = sequence_length
            if sequence_length - start < min_contig:
                # skip last bit of sequence
                break
        else:
            end = start + contig_size
        split_name = '>{}_{}_{}_{}'.format(seq_name, order_number, start, end)
        split_seq = sequence[start:end]

        _ = split_buffer.write('{}\n{}\n'.format(split_name, split_seq))

        order_number += 1
        last_end = end
    assert order_number > 0, 'No split sequences generated: {}'.format(seq_name)
    return split_buffer


rule create_random_mock_assembly:
    input:
        ref_fasta = os.path.join(REFERENCE_FOLDER, 'T2Tv11_T2TC_chm13.fasta')
    output:
        fasta = 'output/assemblies/T2Tv11_randsplit.fasta'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    run:
        with open(output.fasta, 'w') as fasta:
            pass

        chrom_name = None
        chrom_seq = ''
        with open(input.ref_fasta, 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    if chrom_name is None:
                        chrom_name = line.strip().strip('>')
                        continue
                    else:
                        seq_splits = rand_split_sequence(
                            chrom_name,
                            chrom_seq,
                            MIN_GAP_SIZE,
                            MAX_GAP_SIZE,
                            MIN_CONTIG_SIZE
                        )
                        with open(output.fasta, 'a') as fasta:
                            _ = fasta.write(seq_splits.getvalue())
                        chrom_name = line.strip().strip('>')
                        chrom_seq = ''
                        continue
                else:
                    chrom_seq += line.strip()
        seq_splits = rand_split_sequence(
            chrom_name,
            chrom_seq,
            MIN_GAP_SIZE,
            MAX_GAP_SIZE,
            MIN_CONTIG_SIZE
            )
        with open(output.fasta, 'a') as fasta:
            _ = fasta.write(seq_splits.getvalue())
    # END OF RUN BLOCK


def segdup_split_sequence(seq_name, sequence, segdups, min_contig_size):

    import io as io
    import random as rand

    assert not segdups.empty, 'No segdups for chromosome: {}'.format(seg_name)

    sequence_length = len(sequence)

    last_end = 0
    order_number = 0

    split_buffer = io.StringIO()

    # merge overlapping SDs to avoid potentially generating too many small splits
    merge_sd = segdups[['chromStart', 'chromEnd', 'fracMatch']].copy()
    merge_sd['sd_num'] = (merge_sd['chromStart'] > merge_sd['chromEnd'].shift().cummax()).cumsum()
    merge_sd = merge_sd.groupby("sd_num").agg({"chromStart":"min", "chromEnd": "max", "fracMatch": "mean"})

    for sd_num, sd_start, sd_end, sd_id in merge_sd.itertuples(index=False):
        if sd_id < 0.98:
            if rand.random() < 0.5:
                continue
        max_offset = int((sd_end - sd_start) * 0.05)
        offset = rand.randint(0, max_offset)
        start = last_end
        end = sd_start + offset
        if end > sequence_length:
            end = sequence_length
        if end - start < min_contig_size:
            continue
        
        split_name = '>{}_{}_{}_{}'.format(seq_name, order_number, start, end)
        split_seq = sequence[start:end]

        _ = split_buffer.write('{}\n{}\n'.format(split_name, split_seq))

        order_number += 1
        last_end = sd_end - offset
        if last_end + min_contig_size > sequence_length:
            break
    assert order_number > 0, 'No split sequences generated: {}'.format(seq_name)
    return split_buffer


rule create_sdplit_mock_assembly:
    input:
        ref_fasta = os.path.join(REFERENCE_FOLDER, 'T2Tv11_T2TC_chm13.fasta'),
        segdups = os.path.join(REFERENCE_FOLDER, 't2t_segdups.tsv.gz')
    output:
        fasta = 'output/assemblies/T2Tv11_sdsplit.fasta'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    run:
        import pandas as pd

        with open(output.fasta, 'w') as fasta:
            pass

        sd = pd.read_csv(input.segdups, sep='\t', header=0)
        sd = sd[['#chrom', 'chromStart', 'chromEnd', 'fracMatch']]
        sd = sd.loc[sd['fracMatch'].round(2) > 0.95, :].copy()
        sd.sort_values(['#chrom', 'chromStart'], inplace=True, ascending=True)

        chrom_name = None
        chrom_seq = ''
        with open(input.ref_fasta, 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    if chrom_name is None:
                        chrom_name = line.strip().strip('>')
                        continue
                    else:
                        seq_splits = segdup_split_sequence(
                            chrom_name,
                            chrom_seq,
                            sd.loc[sd['#chrom'] == chrom_name, :],
                            MIN_CONTIG_SIZE
                        )
                        with open(output.fasta, 'a') as fasta:
                            _ = fasta.write(seq_splits.getvalue())
                        chrom_name = line.strip().strip('>')
                        chrom_seq = ''
                        continue
                else:
                    chrom_seq += line.strip()
        seq_splits = rand_split_sequence(
            chrom_name,
            chrom_seq,
            sd.loc[sd['#chrom'] == chrom_name, :],
            MIN_CONTIG_SIZE
        )
        with open(output.fasta, 'a') as fasta:
            _ = fasta.write(seq_splits.getvalue())
    # END OF RUN BLOCK


rule count_assembly_kmers:
    input:
        fasta = 'output/assemblies/{assembly}.fasta'
    output:
        kmer_db = directory('output/kmer_db/{assembly}.k15.no-hpc.db'),
        rep_kmer = 'output/kmer_db/{assembly}.k15.no-hpc.rep-grt09998.txt.gz'
    benchmark:
        'rsrc/output/kmer_db/{assembly}.k15.no-hpc.count-dump.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        mem_total_gb = lambda wildcards, attempt: 32 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'meryl count k=15 threads={threads} memory={resources.mem_total_gb} output {output.kmer_db} {input.fasta} && '
        'meryl print greater-than distinct=0.9998 {output.kmer_db} | pigz -p {threads} --best > {output.rep_kmer}'


rule compute_fasta_index:
    input:
        '{folder}/{filename}.fasta'
    output:
        '{folder}/{filename}.fasta.fai'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'samtools faidx {input}'