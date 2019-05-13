
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
        '{filepath}'
    output:
        '{filepath}.crai'
    threads: 8
    shell:
        "samtools index -@ {threads} {input.cram}"


rule samtools_index_bam_alignment:
    input:
        '{filepath}'
    output:
        '{filepath}.bai'
    threads: 8
    shell:
        "samtools index -@ {threads} {input}"


rule validate_fastq:
    input:
        '{filename}.fastq.gz'
    output:
        'output/fastq_validation/{filename}.stats.pck'
    log: 'log/fastq_validation/{filename}.stats.log'
    wildcard_constraints:
        filename="[A-Za-z0-9_\-\.]+"
    threads: 8
    params:
        scriptdir = config['script_dir']
    run:
        exec = '{params.scriptdir}/collect_read_stats.py --debug'
        exec += ' --fastq-input {input} --output {output}'
        exec += ' --chunk-size 25000 --validate'
        exec += ' --num-cpu {threads} &> {log}'
        shell(exec)