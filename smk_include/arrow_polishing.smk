
onsuccess:
    body_text = "Nothing to report\n{log}"
    if config['notify']:
        shell('mail -s "[Snakemake] DGA-Arrow - success" {} <<< "{}"'.format(config['notify_email'], body_text))

onerror:
    if config['notify']:
        shell('mail -s "[Snakemake] DGA-Arrow - ERRROR" {} < {{log}}'.format(config['notify_email']))

rule arrow_polishing_master:
    input:
        expand('output/assembly_polishing/arrow_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.arrow-pass{round}.ctg.fa',
                read_sample=['HG00514.sequel2.pb-clr-25'],
                contig_sample=['HG00514.sequel2.pb-clr-25'],
                haplotype=['tag-h1-un', 'tag-h2-un'],
                round=[1])


checkpoint create_contig_split_files:
    input:
        contigs = 'output/haplotype_assembly/consensus/{contig_sample}.{haplotype}.ctg.fa',
        contig_index = 'output/haplotype_assembly/consensus/{contig_sample}.{haplotype}.ctg.fa.fai',
    output:
        split_files = directory('output/assembly_polishing/contig_splits/{contig_sample}.{haplotype}')
    wildcard_constraints:
        contig_sample = '[\w\.\-]+',
        haplotype = 'tag\-h[12]\-un'
    run:
        import pysam

        output_dir = output.split_files
        os.makedirs(output_dir, exist_ok=True)

        with pysam.FastaFile(input.contigs) as fasta:
            for ref in fasta.references:
                length = fasta.get_reference_length(ref)

                split_name = '{}@0-{}.txt'.format(ref, length)
                split_path = os.path.join(output_dir, split_name)
                with open(split_path, 'w') as split:
                    pass


rule arrow_contig_polishing_pass1:
    input:
        contigs = 'output/haplotype_assembly/consensus/{contig_sample}.{haplotype}.ctg.fa',
        contig_index = 'output/haplotype_assembly/consensus/{contig_sample}.{haplotype}.ctg.fa.fai',
        read_contig_aln = 'output/assembly_polishing/read_contig_alignments/{read_sample}_to_{contig_sample}.{haplotype}.arrow-pass1.bam',
        read_contig_aln_idx =  'output/assembly_polishing/read_contig_alignments/{read_sample}_to_{contig_sample}.{haplotype}.arrow-pass1.bam.pbi',
        contig_split = 'output/assembly_polishing/contig_splits/{contig_sample}.{haplotype}/{ctg_split}.txt'
    output:
        temp('output/assembly_polishing/contig_splits/{contig_sample}.{haplotype}/{read_sample}_to_{contig_sample}.{haplotype}.arrow_pass1.{ctg_split}.ctg.fa')
    log:
        'log/output/assembly_polishing/contig_splits/{contig_sample}.{haplotype}/{read_sample}_to_{contig_sample}.{haplotype}.arrow_pass1.{ctg_split}.log'
    benchmark:
        'run/output/assembly_polishing/contig_splits/{contig_sample}.{haplotype}/{read_sample}_to_{contig_sample}.{haplotype}.arrow_pass1.{ctg_split}.rsrc'
    conda:
        "../environment/conda/conda_pbtools.yml"
    threads: 6
    wildcard_constraints:
        ctg_split = 'ctg[0-9]+@0\-[0-9]+'
    params:
        ref_window = lambda wildcards, input: os.path.basename(input.contig_split).split('.')[0].replace('@', ':')
    shell:
        'variantCaller --algorithm=arrow --log-file {log} --log-level INFO -j {threads} --algorithm arrow --referenceWindow "{params.ref_window}" --reference {input.contigs} -o {output} {input.read_contig_aln}'


def collect_polished_splits(wildcards):
    checkpoint_output = checkpoints.create_contig_split_files.get(**wildcards).output[0]

    split_files = expand('output/assembly_polishing/contig_splits/{contig_sample}.{haplotype}/{read_sample}_to_{contig_sample}.{haplotype}.arrow_pass{round}.{ctg_split}.ctg.fa',
                        contig_sample=wildcards.contig_sample,
                        read_sample=wildcards.read_sample,
                        haplotype=wildcards.haplotype,
                        round=wildcards.round,
                        ctg_split=glob_wildcards(os.path.join(checkpoint_output, '{ctg_split}.txt')).ctg_split)
    return split_files


rule merge_polished_contig_splits:
    input:
        collect_polished_splits
    output:
        'output/assembly_polishing/arrow_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.arrow-pass{round}.ctg.fa'
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

        for contig in input:
            in_buffer.write(open(contig, 'r').read())

        with open(output[0], 'w') as dump:
            _ = dump.write(in_buffer.getvalue())