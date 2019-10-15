
include: 'canonical_dga.smk'
include: 'strandseq_dga.smk'
include: 'aux_utilities.smk'
include: 'run_alignments.smk'

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
        contigs = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}.{hap}.fasta',
        contig_index = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}.{hap}.fasta.fai',
        alignments = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{reference}.{hap}.arrow-p1.psort.pbn.bam',
        aln_index = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{reference}.{hap}.arrow-p1.psort.pbn.bam.pbi',
        sequence_file = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}.{hap}/sequences/{sequence}.seq'
    output:
        'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/split_by_sequence/{hap_reads}.{hap}.arrow-p1.{sequence}.fasta'
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/split_by_sequence/{hap_reads}.{hap}.arrow-p1.{sequence}.log'
    benchmark:
        'run/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/split_by_sequence/{hap_reads}.{hap}.arrow-p1.{sequence}.rsrc'
    conda:
        config['conda_env_pbtools']
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 98304
    params:
        ref_window = lambda wildcards, input: '{}:0-{}'.format(*(open(input.sequence_file).readline().split()[:2]))
    shell:
        'variantCaller --algorithm=arrow --log-file {log} --log-level INFO -j {threads} ' \
            '--referenceWindow "{params.ref_window}" --reference {input.contigs} -o {output} {input.alignments}'


def collect_polished_splits(wildcards):

    checkpoint_output = checkpoints.create_haploid_assembly_sequence_files.get(**wildcards).output[0]

    split_files = expand('output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/split_by_sequence/{hap_reads}.{hap}.arrow-p1.{sequence}.fasta',
                            var_caller=wildcards.var_caller,
                            gq=wildcards.gq,
                            dp=wildcards.dp,
                            reference=wildcards.reference,
                            vc_reads=wildcards.vc_reads,
                            sts_reads=wildcards.sts_reads,
                            pol_reads=wildcards.pol_reads,
                            hap_reads=wildcards.hap_reads,
                            hap=wildcards.hap,
                            sequence=glob_wildcards(os.path.join(checkpoint_output, '{sequence}.seq')).sequence)
    return sorted(split_files)


rule merge_arrow_polished_sequence_splits:
    input:
        collect_polished_splits
    output:
        'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/{hap_reads}.{hap}.arrow-p1.fasta'
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/{hap_reads}.{hap}.arrow-p1.log'
    resources:
        mem_per_cpu_mb = 4096
    run:
        # note to future me:
        # stuff like that runs into trouble because
        # of character limits of single command lines
        # (something related to this quantity:
        # $ getconf ARG_MAX -> 2097152 )
        # Maybe xargs could be a way to circumvent this,
        # but for now, make it naive and slow...
        # DOES NOT WORK:
        #
        #'cat {input} > {output}'

        import io
        in_buffer = io.StringIO()

        with open(log[0], 'w') as logfile:
            total_splits = len(input)
            _ = logfile.write('Start merging operation for {} splits\n'.format(total_splits))
            for num, split in enumerate(input, start=1):
                _ = logfile.write('{}/{}\n'.format(num, total_splits))
                _ = logfile.write('Processing split {}\n'.format(split))
                in_buffer.write(open(split, 'r').read())
                if num % 100 == 0:
                    logfile.flush()

            _ = logfile.write('Dumping merged sequence to output: {}\n'.format(output[0]))
            with open(output[0], 'w') as dump:
                _ = dump.write(in_buffer.getvalue())