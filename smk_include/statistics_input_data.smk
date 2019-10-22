
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
    threads: config['num_cpu_low']
    params:
        script_dir = config['script_dir']
    shell:
        '{params.script_dir}/collect_read_stats.py --debug --input-files {input} --output {output}' \
            ' --chunk-size 250000 --validate --num-cpu {threads} &> {log}'


rule plot_fastq_input_statistics:
    input:
        'output/statistics/fastq_input/stat_dumps/{sample}.stats.pck'
    output:
        'output/plotting/statistics/fastq_input/{sample}.stats.pdf'
    params:
        script_dir = config['script_dir'],
        lower_bound = 7000
    run:
        upper_bounds = {
            'ccs': 25000,
            'clr': 100000,
            'ul': 100000
            }
        tech = wildcards.sample.split('_')[2].split('-')[-1]
        upper_bound = upper_bounds[tech]

        step_sizes = {
            'ccs': 500,
            'clr': 1000,
            'ul': 1000
            }
        step_size = step_sizes[tech]

        exec = '{params.script_dir}/plot_sample_stats.py'
        exec += ' --pck-input {input} --genome-length 3G'
        exec += ' --sample-name {wildcards.sample}'
        exec += ' --lowest-bin ' + str(params.lower_bound)
        exec += ' --highest-bin ' + str(upper_bound)
        exec += ' --step-size ' + str(step_size)
        exec += ' --output {output}'
        shell(exec)


rule compute_statistics_haplosplit_fastq:
    input:
        'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.fastq.gz'
    output:
        'output/statistics/fastq_haplosplit/stat_dumps/{approach}/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.stats.pck'
    log: 'log/output/statistics/fastq_haplosplit/stat_dumps/{approach}/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.log'
    benchmark: 'run/output/statistics/fastq_haplosplit/stat_dumps/{approach}/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.rsrc'
    threads: config['num_cpu_low']
    params:
        script_dir = config['script_dir']
    shell:
        '{params.script_dir}/collect_read_stats.py --debug --input-files {input} --output {output}' \
            ' --chunk-size 250000 --validate --num-cpu {threads} &> {log}'


rule plot_fastq_haplosplit_statistics:
    input:
        'output/statistics/fastq_haplosplit/stat_dumps/{approach}/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.stats.pck'
    output:
        'output/plotting/statistics/fastq_haplosplit/{approach}/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.stats.pdf'
    params:
        script_dir = config['script_dir'],
        lower_bound = 7000
    run:
        upper_bounds = {
            'ccs': 25000,
            'clr': 100000,
            'ul': 100000
            }
        tech = wildcards.sample.split('_')[2].split('-')[-1]
        upper_bound = upper_bounds[tech]

        step_sizes = {
            'ccs': 500,
            'clr': 1000,
            'ul': 1000
            }
        step_size = step_sizes[tech]

        exec = '{params.script_dir}/plot_sample_stats.py'
        exec += ' --pck-input {input} --genome-length 3G'
        exec += ' --sample-name {wildcards.sample}'
        exec += ' --lowest-bin ' + str(params.lower_bound)
        exec += ' --highest-bin ' + str(upper_bound)
        exec += ' --step-size ' + str(step_size)
        exec += ' --text-size 14'
        exec += ' --output {output}'
        shell(exec)

