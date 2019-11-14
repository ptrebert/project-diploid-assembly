
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
    resources:
        mem_total_mb = 24576,  # dependent on chunk size and avg. read length
        mem_per_cpu_mb = int(24576 / config['num_cpu_low']),
        runtime_hrs=lambda wildcards: 1 if '-ccs' in wildcards.filename else 2
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


rule filter_squashed_assembly_by_size:
    input:
        'references/assemblies/{sample}_sqa-{assembler}.fasta'
    output:
        fasta = 'references/assemblies/filtered/{sample}_sqa-{assembler}-100kb.fasta',
        stats = 'output/statistics/assemblies/{sample}_sqa-{assembler}-100kb.stats.tsv'
    log:
        'log/references/assemblies/{sample}_sqa-{assembler}-100kb.log'
    params:
        scriptdir = config['script_dir'],
        min_contig_size = config['min_contig_size'],
    shell:
        '{params.scriptdir}/filter_squashed_assembly.py --debug --input-fasta {input}' \
            ' --output-fasta {output.fasta} --output-metrics {output.stats}' \
            ' --min-size {params.min_contig_size} &> {log}'


rule compute_statistics_haplosplit_fastq:
    input:
        'output/diploid_assembly/strandseq_joint/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.fastq.gz'
    output:
        'output/statistics/fastq_haplosplit/stat_dumps/strandseq_joint/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.stats.pck'
    log:
        'log/output/statistics/fastq_haplosplit/stat_dumps/strandseq_joint/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.log'
    benchmark:
        'run/output/statistics/fastq_haplosplit/stat_dumps/strandseq_joint/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.rsrc'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = 24576,  # dependent on chunk size and avg. read length
        mem_per_cpu_mb = int(24576 / config['num_cpu_low']),
        runtime_hrs=lambda wildcards: 1 if '-ccs' in wildcards.filename else 2
    params:
        script_dir = config['script_dir']
    shell:
        '{params.script_dir}/collect_read_stats.py --debug --input-files {input} --output {output}' \
            ' --chunk-size 250000 --validate --num-cpu {threads} &> {log}'


rule plot_fastq_haplosplit_statistics:
    input:
        'output/statistics/fastq_haplosplit/stat_dumps/strandseq_joint/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.stats.pck'
    output:
        'output/plotting/statistics/fastq_haplosplit/strandseq_joint/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.stats.pdf'
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

