
include: 'preprocess_input.smk'

localrules: master_statistics_input_data

rule master_statistics_input_data:
    input:


rule validate_complete_input_fastq:
    input:
        'input/fastq/complete/{filename}.fastq.gz'
    output:
        'output/statistics/fastq_input/stat_dumps/{filename}.stats.pck'
    log: 'log/output/statistics/fastq_input/stat_dumps/{filename}.stats.log'
    benchmark: 'run/output/statistics/fastq_input/stat_dumps/{filename}.stats.rsrc'
    threads: 4
    params:
        scriptdir = config['script_dir']
    run:
        exec = '{params.scriptdir}/collect_read_stats.py --debug'
        exec += ' --input-files {input} --output {output}'
        exec += ' --chunk-size 250000 --validate'
        exec += ' --num-cpu {threads} &> {log}'
        shell(exec)