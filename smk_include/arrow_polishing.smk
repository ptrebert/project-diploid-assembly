
include: 'canonical_dga.smk'
include: 'strandseq_dga_joint.smk'
include: 'strandseq_dga_split.smk'
include: 'aux_utilities.smk'
include: 'run_alignments.smk'


rule arrow_contig_polishing_pass1:
    input:
        contigs = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}.{sequence}.fasta',
        seq_info = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fasta/{hap_reads}-{assembler}.{hap}.{sequence}.fasta.fai',
        alignments = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.psort.pbn.bam',
        aln_index = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/alignments/{pol_reads}_map-to_{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.psort.pbn.bam.pbi',
    output:
        'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.fasta'
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.log'
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.rsrc'
    conda:
        config['conda_env_pbtools']
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = int(32768 / config['num_cpu_medium']),
        mem_total_mb = 32768,
        runtime_hrs = 12
    shell:
        'variantCaller --algorithm=arrow --log-file {log} --log-level INFO -j {threads} ' \
            ' --reference {input.contigs} -o {output} {input.alignments}'


def collect_arrow_polished_contigs(wildcards):

    reference_folder = os.path.join('output/reference_assembly/clustered', wildcards.sts_reads)
    seq_output_dir = checkpoints.create_assembly_sequence_files.get(folder_path=reference_folder,
                                                                    reference=wildcards.reference).output[0]

    checkpoint_wildcards = glob_wildcards(
        os.path.join(seq_output_dir, '{sequence}.seq')
        )

    polished_contigs = expand(
        'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.{sequence}.arrow-p1.fasta',
        var_caller=wildcards.var_caller,
        qual=wildcards.qual,
        gq=wildcards.gq,
        reference=wildcards.reference,
        vc_reads=wildcards.vc_reads,
        sts_reads=wildcards.sts_reads,
        pol_reads=wildcards.pol_reads,
        hap_reads=wildcards.hap_reads,
        assembler=wildcards.assembler,
        hap=wildcards.hap,
        sequence=checkpoint_wildcards.sequence
    )

    return polished_contigs


rule write_arrow_polished_contigs_fofn:
    input:
        contigs = collect_arrow_polished_contigs
    output:
        fofn = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.arrow-p1.fofn'
    resources:
        runtime_hrs = 0,
        runtime_min = 10
    run:
        contigs = sorted(collect_arrow_polished_contigs(wildcards))

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(contigs):
                if not os.path.isfile(file_path):
                    if os.path.isdir(file_path):
                        # this is definitely wrong
                        raise AssertionError('Expected file path for Arrow polished FASTA, but received directory: {}'.format(file_path))
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule merge_arrow_polished_contigs:
    input:
        fofn = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.arrow-p1.fofn'
    output:
        fofn = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_fasta/{hap_reads}-{assembler}.{hap}.arrow-p1.fasta'
    resources:
        runtime_hrs = 0,
        runtime_min = 20
    params:
        polished_fastas = lambda wildcards, input: load_fofn_file(input)
    shell:
        'cat {params.polished_fastas} > {output}'