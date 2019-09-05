
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