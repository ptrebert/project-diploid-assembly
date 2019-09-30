
rule compute_md5_checksum:
    input:
        '{filepath}'
    output:
        '{filepath}.md5'
    shell:
        "md5sum {input} > {output}"


rule samtools_index_cram_alignment:
    input:
        cram = '{filepath}'
    output:
        crai = '{filepath}.crai'
    threads: 8
    shell:
        "samtools index -@ {threads} {input.cram}"


rule samtools_index_bam_alignment:
    input:
        bam = '{filepath}'
    output:
        bai = '{filepath}.bai'
    threads: 8
    shell:
        "samtools index -@ {threads} {input.bam}"


rule samtools_convert_sam_to_bam:
    input:
        sam = '{filepath}.sam'
    output:
        bam = '{filepath}.sam.bam'
    wildcard_constraints:
        filepath = '[\w\-\/]+'
    threads: 8
    shell:
        "samtools view -o {output.bam} -b -@ {threads} {input.sam}"


rule samtools_position_sort_bam_alignment:
    input:
        unsorted_bam = '{filepath}.sam.bam'
    output:
        sorted_bam = '{filepath}.psort.sam.bam'
    wildcard_constraints:
        filepath = '[\w\-\/]+'
    threads: 8
    resources:
        mem_mb = 8 * (20 * 1024)  # eight times twenty gigabyte
    params:
        mem_per_thread = '20G'
    run:
        exec = 'samtools sort'
        exec += ' -m {params.mem_per_thread}'
        exec += ' --threads {threads}'
        exec += ' -o {output.sorted_bam}'
        exec += ' {input.unsorted_bam}'
        shell(exec)


rule samtools_index_fasta:
    input:
        '{filepath}.fasta'
    output:
        '{filepath}.fasta.fai'
    shell:
        'samtools faidx {input}'


rule bgzip_file_copy:
    input:
        '{filepath}'
    output:
        '{filepath}.bgz'
    shell:
        'bgzip -c {input} > {output}'


rule bcftools_index_bgzipped_file:
    input:
        '{filepath}.bgz'
    output:
        '{filepath}.bgz.tbi'
    shell:
        'bcftools index --tbi {input}'


checkpoint create_assembly_sequence_files:
    input:
        'references/assemblies/{reference}.fasta.fai'
    output:
        directory('references/assemblies/{reference}/sequences')
    wildcard_constraints:
        reference = '[\-\w]+'
    run:
        output_dir = output[0]
        os.makedirs(output_dir, exist_ok=True)
        with open(input[0], 'r') as fai:
            for line in fai:
                seq_name = line.split('\t')[0]
                output_path = os.path.join(output_dir, seq_name + '.seq')
                with open(output_path, 'w') as dump:
                    _ = dump.write(line)


def collect_assembly_sequence_files(wildcards):
    """
    """
    seq_output_dir = checkpoints.create_assembly_sequence_files.get(reference=wildcards.reference).output[0]

    checkpoint_wildcards = glob_wildcards(seq_output_dir, '{sequence}.seq')

    return expand(os.path.join(seq_output_dir, '{sequence}.seq'),
                    sequence=checkpoint_wildcards.sequence)


rule generate_bwa_index:
    input:
        reference = 'references/assemblies/{reference}.fasta'
    output:
        'references/assemblies/bwa_index/{reference}.amb',
        'references/assemblies/bwa_index/{reference}.ann',
        'references/assemblies/bwa_index/{reference}.bwt',
        'references/assemblies/bwa_index/{reference}.pac',
        'references/assemblies/bwa_index/{reference}.sa'
    params:
        prefix = lambda wildcards, output: output[0].split('.')[0]
    log:
        'log/references/assemblies/bwa_index/{reference}.log'
    benchmark:
        'run/references/assemblies/bwa_index/{reference}.rsrc'
    shell:
        'bwa index -p {params.prefix} {input.reference} &> {log}'
