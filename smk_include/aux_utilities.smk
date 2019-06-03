
rule compute_md5_checksum:
    input:
        '{filepath}'
    output:
        '{filepath}.md5'
    shell:
        "md5sum {input} > {output}"


#rule gzip_file_copy:
#    input:
#        '{filepath}'
#    output:
#        '{filepath}.gz'
#    shell:
#        "gzip -c {input} > {output}"


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
        exec += ' --fastq-input {input} --output {output}'
        exec += ' --chunk-size 25000 --validate'
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
        exec += ' --fastq-input {input} --output {output}'
        exec += ' --chunk-size 25000 --validate'
        exec += ' --num-cpu {threads} &> {log}'
        shell(exec)