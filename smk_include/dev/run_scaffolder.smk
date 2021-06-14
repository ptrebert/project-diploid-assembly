
REFERENCE_FOLDER = '/beeond/data/references'

MIN_CONTIG_SIZE = 20000
MIN_GAP_SIZE = 1000
MAX_GAP_SIZE = 500000

# filter: select informative alignments


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


rule master:
    input:
        'output/assemblies/T2Tv11_randsplit.fasta.fai',
        'output/assemblies/T2Tv11_sdsplit.fasta.fai',
        'output/graphs/ont_to_assm/chm13_ONTrel7_MAP-TO_T2Tv11_randsplit.wmap-k15.filt.noseq.gfa',
        'output/graphs/ont_to_assm/chm13_ONTrel7_MAP-TO_T2Tv11_sdsplit.wmap-k15.filt.noseq.gfa',

def rand_split_sequence(seq_name, sequence, min_gap, max_gap, min_contig):

    import io
    import random as rand

    if seq_name in ['chrM', 'chrY']:
        return io.StringIO()

    sequence_length = len(sequence)

    max_contig = sequence_length // 4
    split_buffer = io.StringIO()

    last_end = 0
    order_number = 0
    while 1:
        if last_end != 0:
            start = last_end + rand.randint(min_gap, max_gap)
        else:
            start = last_end
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
        ref_fasta = ancient(os.path.join(REFERENCE_FOLDER, 'T2Tv11_T2TC_chm13.fasta'))
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

    if seq_name in ['chrM', 'chrY']:
        return io.StringIO()

    assert not segdups.empty, 'No segdups for chromosome: {}'.format(seg_name)

    sequence_length = len(sequence)

    last_end = 0
    order_number = 0

    split_buffer = io.StringIO()

    # merge overlapping SDs to avoid potentially generating too many small splits
    merge_sd = segdups[['chromStart', 'chromEnd', 'fracMatch']].copy()
    merge_sd['sd_num'] = (merge_sd['chromStart'] > merge_sd['chromEnd'].shift().cummax()).cumsum()
    merge_sd = merge_sd.groupby("sd_num").agg({"chromStart":"min", "chromEnd": "max", "fracMatch": "mean"})

    for sd_start, sd_end, sd_id in merge_sd.itertuples(index=False):
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
        ref_fasta = ancient(os.path.join(REFERENCE_FOLDER, 'T2Tv11_T2TC_chm13.fasta')),
        segdups = ancient(os.path.join(REFERENCE_FOLDER, 't2t_segdups.tsv.gz'))
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
        seq_splits = segdup_split_sequence(
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
        rep_kmer = 'output/kmer_db/{assembly}.k15.no-hpc.rep-grt09998.txt'
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
        'meryl print greater-than distinct=0.9998 {output.kmer_db} > {output.rep_kmer}'


rule compute_fasta_index:
    input:
        '{folder}/{filename}.fasta'
    output:
        '{folder}/{filename}.fasta.fai'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'samtools faidx {input}'


rule align_ont_reads_to_assembly:
    input:
        reads = '/beeond/data/ont/chm13/rel7.fastq.gz',
        reference = 'output/assemblies/{assembly}.fasta',
        ref_repkmer = 'output/kmer_db/{assembly}.k15.no-hpc.rep-grt09998.txt',
    output:
        paf = 'output/alignments/ont_to_assm/chm13_ONTrel7_MAP-TO_{assembly}.wmap-k15.paf.gz'
    log:
        'log/output/alignments/ont_to_assm/chm13_ONTrel7_MAP-TO_{assembly}.wmap-k15.log'
    benchmark:
        'rsrc/output/alignments/ont_to_assm/chm13_ONTrel7_MAP-TO_{assembly}.wmap-k15.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_max']
    resources:
        mem_total_mb = lambda wildcards, attempt: 131072 * attempt,
        runtime_hrs = lambda wildcards, attempt: 167 * attempt,
    params:
        align_threads = config['num_cpu_max'] - 2,
        zip_threads = 2
    shell:
        'winnowmap -W {input.ref_repkmer} -k 15 -t {params.align_threads} -x map-ont '
        '--secondary=no {input.reference} {input.reads} 2> {log} '
        '|'
        ' pigz -p {params.zip_threads} --best > {output.paf}'


rule extract_usable_alignments:
    input:
        paf = 'output/alignments/ont_to_assm/chm13_ONTrel7_MAP-TO_{assembly}.wmap-k15.paf.gz'
    output:
        paf = 'output/alignments/ont_to_assm/chm13_ONTrel7_MAP-TO_{assembly}.wmap-k15.filt.paf'
    conda:
        '../../environment/conda/conda_pyscript.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 6144 * attempt
    params:
        script_exec = lambda wildcards: find_script_path('scf_filter_paf.py', 'scaffolding'),
        mmr = 300,
        mmq = 60,
        tip = 5,
        cont = 90,
    shell:
        '{params.script_exec} --paf {input.paf} --output {output.paf} '
            '--min-matched-residues {params.mmr} --min-mapq-threshold {params.mmq} '
            '--tip-boundary {params.tip} --containment {params.cont}'


def determine_read_contig_order(record):
    """
    Since the overlap-relation between segments
    implies an order, need to determine who comes
    first in a Link line
    """
    contig_end = 'left'
    # check simply to which end the alignment start is closer
    dist_left = int(record['aln_tstart'])
    dist_right = int(record['target_length']) - int(record['aln_tend'])
    if dist_right < dist_left:
        contig_end = 'right'
    return contig_end    


rule convert_paf_to_gfa:
    """
    this is a temporary solution to visualize the graph/alignment underlying
    the scaffolding (the input). This converts to GFAv1 mainly for reasons
    of compatibility with Bandage. GFAv2 would be a better fit to precisely record
    the information in the PAF file but has lower tool support.
    To be re-done in form of a script if visualization issue is solved.
    """
    input:
        paf = 'output/alignments/ont_to_assm/chm13_ONTrel7_MAP-TO_{assembly}.wmap-k15.filt.paf'
    output:
        gfa = 'output/graphs/ont_to_assm/chm13_ONTrel7_MAP-TO_{assembly}.wmap-k15.filt.noseq.gfa'
    run:
        import io
        import csv

        gfa_segments = io.StringIO()
        segment_line = 'S\t{name}\t*\tLN:i:{length}\n'

        gfa_links = io.StringIO()
        link_line = 'L\t{seg1}\t{orient1}\t{seg2}\t{orient1}\t{overlap}M\tNM:i:{mismatch}\tID:Z:{ln}\n'

        segments_added = set()

        with open(input.paf, 'r', newline='') as paf:
            columns = paf.readline().strip().split()
            paf_records = csv.DictReader(paf, fieldnames=columns, delimiter='\t')
            for ln, record in enumerate(paf_records, start=2):
                read_segment = record['query_name'].replace('-', '_')  # cannot use - in names
                if read_segment not in segments_added:
                    gfa_segments.write(segment_line.format(**{'name': read_segment, 'length': record['query_length']}))
                    segments_added.add(read_segment)
                if record['target_name'] not in segments_added:
                    gfa_segments.write(segment_line.format(**{'name': record['target_name'], 'length': record['target_length']}))
                    segments_added.add(record['target_name'])
                # alignments are unique
                contig_end = determine_read_contig_order(record)
                if contig_end == 'left':
                    # read before contig
                    seg1 = read_segment
                    orient1 = record['aln_qstrand']
                    seg2 = record['target_name']
                    orient2 = '+'
                else:
                    seg2 = read_segment
                    orient2 = record['aln_qstrand']
                    seg1 = record['target_name']
                    orient1 = '+'
                # NM: number of mismatches/gaps
                # estimate by subtracting matches
                # from alignment block length
                mismatches = int(record['aln_block_length']) - int(record['aln_num_match'])
                gfa_links.write(link_line.format(**{
                    'seg1': seg1,
                    'orient1': orient1,
                    'seg2': seg2,
                    'orient2': orient2,
                    'overlap': record['aln_num_match'],  # Bandage does not like "*" CIGAR strings
                    'mismatch': mismatches,
                    'ln': ln
                }))
        
        with open(output.gfa, 'w') as gfa:
            _ = gfa.write('H\tVN:Z:GFAv1\n')
            _ = gfa.write(gfa_segments.getvalue())
            _ = gfa.write(gfa_links.getvalue())
    # END OF RUN BLOCK

    