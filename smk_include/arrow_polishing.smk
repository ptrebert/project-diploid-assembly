

checkpoint create_haploid_assembly_sequence_files:
    input:
        'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}.{hap}.fasta.fai'
    output:
        directory('output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}.{hap}/sequences')
    run:
        output_dir = output[0]
        os.makedirs(output_dir, exist_ok=True)
        with open(input[0], 'r') as fai:
            for line in fai:
                seq_name = line.split('\t')[0]
                output_path = os.path.join(output_dir, seq_name + '.seq')
                with open(output_path, 'w') as dump:
                    _ = dump.write(line)


rule arrow_contig_polishing_pass1:
    input:
        contigs = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}.{hap}.fasta,
        contig_index = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}.{hap}.fasta.fai',
        alignments = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{reference}.{hap}.arrow-p1.psort.pbn.bam',
        aln_index = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{reference}.{hap}.arrow-p1.psort.pbn.bam.bai',
        sequence_file = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}.{hap}/sequences/{sequence}.seq'
    output:
        'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/{hap_reads}.{hap}.arrow-p1.{sequence}.fasta'
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/{hap_reads}.{hap}.arrow-p1.{sequence}.log'
    benchmark:
        'run/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/{hap_reads}.{hap}.arrow-p1.{sequence}.rsrc'
    conda:
        "../environment/conda/conda_pbtools.yml"
    threads: 8
    params:
        ref_window = lambda wildcards, input: '{}:0-{}'.format(open(input.sequence_file).readline().split()[:2])
    shell:
        'variantCaller --algorithm=arrow --log-file {log} --log-level INFO -j {threads} ' \
            '--referenceWindow "{params.ref_window}" --reference {input.contigs} -o {output} {input.alignments}'


#def collect_polished_splits(wildcards):
#    checkpoint_output = checkpoints.create_contig_split_files.get(**wildcards).output[0]
#
#    split_files = expand('output/assembly_polishing/contig_splits/{contig_sample}.{haplotype}/{read_sample}_to_{contig_sample}.{haplotype}.arrow_pass{round}.{ctg_split}.ctg.fa',
#                        contig_sample=wildcards.contig_sample,
#                        read_sample=wildcards.read_sample,
#                        haplotype=wildcards.haplotype,
#                        round=wildcards.round,
#                        ctg_split=glob_wildcards(os.path.join(checkpoint_output, '{ctg_split}.txt')).ctg_split)
#    return split_files
#
#
#rule merge_polished_contig_splits:
#    input:
#        collect_polished_splits
#    output:
#        'output/assembly_polishing/arrow_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.arrow-pass{round}.ctg.fa'
#    run:
#        # note to future me:
#        # stuff like that runs into trouble because
#        # of character limits of single command lines
#        # (something related to this quantity:
#        # $ getconf ARG_MAX -> 2097152 )
#        # Maybe xargs could be a way to circumvent this,
#        # but for now, make it naive and slow...
#        # DOES NOT WORK:
#        #
#        #'cat {input} > {output}'
#
#        import io
#        in_buffer = io.StringIO()
#
#        for contig in input:
#            in_buffer.write(open(contig, 'r').read())
#
#        with open(output[0], 'w') as dump:
#            _ = dump.write(in_buffer.getvalue())