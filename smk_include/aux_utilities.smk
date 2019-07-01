
rule compute_md5_checksum:
    input:
        '{filepath}'
    output:
        '{filepath}.md5'
    shell:
        "md5sum {input} > {output}"


rule gzip_file_copy:
    input:
        '{filepath}'
    output:
        '{filepath}.gz'
    shell:
        "gzip -c {input} > {output}"


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


rule samtools_sort_bam_alignment:
    input:
        unsorted_bam = '{filepath}.bam'
    output:
        sorted_bam = '{filepath}.sort.bam'
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
        '{filepath}/{filename}.fa'
    output:
        '{filepath}/{filename}.fa.fai'
    wildcard_constraints:
        filename = '[A-Za-z0-9\.\-_]+'
    shell:
        'samtools faidx {input}'


rule validate_split_fastq:
    input:
        'output/haplotype_partitioning/splitting/{filename}.fastq.gz'
    output:
        'output/fastq_validation/{filename}.stats.pck'
    log: 'log/fastq_validation/{filename}.stats.log'
    benchmark: 'run/fastq_validation/{filename}.stats.rsrc'
    threads: 4
    params:
        scriptdir = config['script_dir']
    run:
        exec = '{params.scriptdir}/collect_read_stats.py --debug'
        exec += ' --input-files {input} --output {output}'
        exec += ' --chunk-size 150000 --validate'
        exec += ' --num-cpu {threads} &> {log}'
        shell(exec)


rule validate_input_fastq:
    input:
        'input/read_data/diploid_assembly_input/{filename}.fastq.gz'
    output:
        'output/fastq_validation/{filename}.stats.pck'
    log: 'log/fastq_validation/{filename}.stats.log'
    benchmark: 'run/fastq_validation/{filename}.stats.rsrc'
    threads: 4
    params:
        scriptdir = config['script_dir']
    run:
        exec = '{params.scriptdir}/collect_read_stats.py --debug'
        exec += ' --input-files {input} --output {output}'
        exec += ' --chunk-size 150000 --validate'
        exec += ' --num-cpu {threads} &> {log}'
        shell(exec)