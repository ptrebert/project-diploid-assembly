
sts_path = '/scratch/bioinf/projects/diploid-genome-assembly/pebert/test_ccs/run_folder/input/fastq/NA12878_eriba_il25k-100pe_sseq'

sts_mate1 = sorted([os.path.join(sts_path, f) for f in os.listdir(sts_path) if f.endswith('_1.fastq.gz')])
sts_mate2 = sorted([os.path.join(sts_path, f) for f in os.listdir(sts_path) if f.endswith('_2.fastq.gz')])

library_ids = sorted([os.path.basename(f).rsplit('_', 1)[0] for f in sts_mate1])

rule master:
    input:
        'output/alignments/NA12878_reads_to_ref.filt.sam.bam',
        'output/readset/NA12878_demo_pbsq2-ccs_1000.fastq.gz',
         expand('output/sts_fastq/{library_id}_{num}.fastq.gz',
                library_id=[l.replace('eriba', 'demo') for l in library_ids],
                num=[1, 2])



rule select_demo_contigs:
    input:
        expand('/scratch/bioinf/users/pebert/data_sources/na12878_ccs_ref_clusters/cluster{num}.fasta',
               num=list(range(1, 24)))
    output:
        'output/reference/NA12878_demo_reference.fasta'
    run:
        import os
        import io
        import collections as col
        split_seq = 'N' * 100
        length_limit = 100 * 1e6
        lower_bound = 5 * 1e6
        upper_bound = 10 * 1e6

        input_cluster = sorted([(os.stat(fp).st_size, fp) for fp in list(input)], reverse=True)

        selected_contig_length = 0
        cluster_contig_counts = col.Counter()
        buffer = io.StringIO()
        for _, fp in input_cluster:
            cluster_select = 0
            with open(fp, 'r') as fasta:
                cluster_name = fasta.readline().strip().strip('>')
                sequence = fasta.read().strip().replace('\n', '')
                contigs = sequence.split(split_seq)
                assert len(contigs) > 1, 'Split sequence is wrong: {}'.format(fp)
                for number, contig_seq in enumerate(contigs, start=1):
                    contig_length = len(contig_seq)
                    if lower_bound < contig_length < upper_bound:
                        if selected_contig_length + contig_length > length_limit:
                            continue
                        seq_header = cluster_name + '_{}_demo_{}'.format(number, contig_length)
                        _ = buffer.write('>' + seq_header + '\n')
                        _ = buffer.write(contig_seq + '\n')
                        cluster_select += 1
                        selected_contig_length += contig_length
                        cluster_contig_counts[cluster_name] += 1
                        if cluster_select == 2:
                            break

        #print(cluster_contig_counts)
        # Counter({'cluster8': 2, 'cluster15': 2, 'cluster5': 2, 'cluster2': 2, 'cluster17': 2, 'cluster11': 2, 'cluster21': 1, 'cluster14': 1, 'cluster12': 1})

        with open(output[0], 'w') as dump:
            _ = dump.write(buffer.getvalue())


rule align_reads_to_demo_reference:
    input:
        ref = 'output/reference/NA12878_demo_reference.fasta',
        reads = '/scratch/bioinf/projects/diploid-genome-assembly/pebert/test_ccs/run_folder/input/fastq/NA12878_giab_pbsq2-ccs_1000.fastq.gz'
    output:
        bam = 'output/alignments/NA12878_reads_to_ref.sam.bam'
    conda:
        '../../environment/conda/conda_pbtools.yml'
    threads: 48
    params:
        tempdir = lambda wildcards: os.path.join('temp', 'pbmm2', 'reads_to_ref')
    shell:
         'TMPDIR={params.tempdir} '
         'pbmm2 align --log-level INFO --sort --sort-memory 16384M --no-bai '
         ' --alignment-threads 40 --sort-threads 8 --preset CCS --min-length 5000 '
         ' --rg "@RG\\tID:NA12878_demo\\tSM:NA12878" '
         ' --sample NA12878 {input.ref} {input.reads} {output.bam}'


rule filter_reads_to_ref_alignment:
    input:
        'output/alignments/NA12878_reads_to_ref.sam.bam'
    output:
        'output/alignments/NA12878_reads_to_ref.filt.sam.bam'
    threads: 8
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'samtools view -b -q 30 -F 3844 -@ 8 -o {output} {input}'


rule dump_reads_to_ref_alignment:
    input:
         'output/alignments/NA12878_reads_to_ref.sam.bam'
    output:
          'output/readset/NA12878_demo_pbsq2-ccs_1000.fastq.gz'
    threads: 8
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         'samtools fastq -@ 8 -n -0 /dev/stdout {input} | gzip > {output}'


rule generate_bwa_index:
    input:
         'output/reference/NA12878_demo_reference.fasta'
    output:
          'output/reference/NA12878_demo_reference/bwa_index/NA12878_demo_reference.bwt'
    conda:
          '../../environment/conda/conda_biotools.yml'
    params:
        prefix = lambda wildcards, output: output[0].rsplit('.', 1)[0]
    shell:
        'bwa index -p {params.prefix} {input}'


rule strandseq_alignments:
    input:
         index = 'output/reference/NA12878_demo_reference/bwa_index/NA12878_demo_reference.bwt',
         mate1 = os.path.join(sts_path, '{library_id}_1.fastq.gz'),
         mate2 = os.path.join(sts_path, '{library_id}_2.fastq.gz')
    output:
          'output/sts_align/{library_id}.bam'
    log:
       bwa = 'log/output/sts_align/{library_id}.bwa.log',
       samtools = 'log/output/sts_align/{library_id}.samtools.log'
    threads: 4
    conda:
         '../../environment/conda/conda_biotools.yml'
    params:
          idx_prefix = lambda wildcards, input: input.index.rsplit('.', 1)[0]
    shell:
         'bwa mem -t {threads}'
         ' -R "@RG\\tID:{wildcards.library_id}\\tPL:Illumina\\tSM:NA12878_demo"'
         ' -v 2 {params.idx_prefix} {input.mate1} {input.mate2} 2> {log.bwa} | '
         ' samtools view -b -F 2308 -q 10 /dev/stdin > {output} 2> {log.samtools}'


rule dump_bam_to_fastq:
    input:
         'output/sts_align/{library_id}.bam'
    output:
          'output/sts_fastq/{library_id}_1.fastq',
          'output/sts_fastq/{library_id}_2.fastq'
    threads: 4
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         'samtools fastq -s /dev/null -N -1 {output[0]} -2 {output[1]} {input}'


rule gzip:
    input:
         'output/sts_fastq/{individual}_eriba_{library_id}.fastq'
    output:
          'output/sts_fastq/{individual}_demo_{library_id}.fastq.gz'
    shell:
         'gzip -c {input} > {output}'
