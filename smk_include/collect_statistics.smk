
localrules: master_statistics_input_data

rule master_statistics_input_data:
    input:
        []


rule compute_statistics_complete_input_fastq:
    input:
        fastq = 'input/fastq/{sample}.fastq.gz',
        faidx = 'references/assemblies/' + config['use_genome_size'] +'.fasta.fai',
    output:
        dump = 'output/statistics/stat_dumps/{sample}.fastq.pck',
        summary = 'input/fastq/{sample}.stats',
    log: 'log/output/statistics/stat_dumps/{sample}.fastq.log',
    benchmark: 'run/output/statistics/stat_dumps/{sample}.fastq.t2.rsrc'
    threads: 2
    resources:
        runtime_hrs= 8,
        mem_total_mb = 4096,
        mem_per_cpu_mb = 2048,
    conda:
         '../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
        '{params.script_exec} --debug --input-files {input.fastq} '
        '--output {output.dump} --summary-output {output.summary} '
        '--copy-stats-dump '
        'output/statistics/stat_dumps/{wildcards.sample}.pbn.bam.pck '
        'output/statistics/stat_dumps/{wildcards.sample}.fasta.pck '
        '--copy-summary '
        'input/bam/{wildcards.sample}.stats '
        'input/fasta/{wildcards.sample}.stats '
        '--num-cpu {threads} --genome-size-file {input.faidx} &> {log}'


rule compute_statistics_complete_input_fasta:
    input:
        fasta = 'input/fasta/{sample}.fasta',
        faidx = 'references/assemblies/' + config['use_genome_size'] +'.fasta.fai',
    output:
        dump = 'output/statistics/stat_dumps/{sample}.fasta.pck',
        summary = 'input/fasta/{sample}.stats'
    log: 'log/output/statistics/stat_dumps/{sample}.fasta.log',
    benchmark: 'run/output/statistics/stat_dumps/{sample}.fasta.t2.rsrc'
    threads: 2
    resources:
        runtime_hrs= 8,
        mem_total_mb = 4096,
        mem_per_cpu_mb = 2048,
    conda:
         '../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
        '{params.script_exec} --debug --input-files {input.fasta} '
        '--output {output.dump} --summary-output {output.summary} '
        '--copy-stats-dump '
        'output/statistics/stat_dumps/{wildcards.sample}.pbn.bam.pck '
        'output/statistics/stat_dumps/{wildcards.sample}.fastq.pck '
        '--copy-summary '
        'input/bam/{wildcards.sample}.stats '
        'input/fastq/{wildcards.sample}.stats '
        '--num-cpu {threads} --genome-size-file {input.faidx} &> {log}'


rule compute_statistics_split_cluster_fastq:
    input:
         fasta = 'output/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq.gz',
         faidx = 'output/reference_assembly/clustered/{sts_reads}/{reference}/sequences/{sequence}.seq',
    output:
          dump = 'output/statistics/stat_dumps/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq.pck',
          summary = 'output/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.stats',
    log: 'log/output/statistics/stat_dumps/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq.log',
    benchmark: 'run/output/statistics/stat_dumps/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/draft/haploid_fastq/{hap_reads}.{hap}.{sequence}.fastq.t2.rsrc'
    threads: 2
    resources:
             runtime_hrs= 1,
             mem_total_mb = lambda wildcards, attempt: 1024 + 1024 * attempt,
             mem_per_cpu_mb = lambda wildcards, attempt: (1024 + 1024 * attempt) // 2,
    conda:
         '../environment/conda/conda_pyscript.yml'
    params:
          script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
         '{params.script_exec} --debug --input-files {input.fasta} '
         '--output {output.dump} --summary-output {output.summary} '
         '--num-cpu {threads} --genome-size-file {input.faidx} &> {log}'


rule compute_statistics_split_cluster_fasta:
    input:
         fasta = 'output/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/draft/haploid_fasta/{hap_reads}.{hap}.{sequence}.fasta',
         faidx = 'output/reference_assembly/clustered/{sts_reads}/{reference}/sequences/{sequence}.seq',
    output:
          dump = 'output/statistics/stat_dumps/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/draft/haploid_fasta/{hap_reads}.{hap}.{sequence}.fasta.pck',
          summary = 'output/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/draft/haploid_fasta/{hap_reads}.{hap}.{sequence}.stats',
    log: 'log/output/statistics/stat_dumps/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/draft/haploid_fasta/{hap_reads}.{hap}.{sequence}.fasta.log',
    benchmark: 'run/output/statistics/stat_dumps/diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/draft/haploid_fasta/{hap_reads}.{hap}.{sequence}.fasta.t2.rsrc'
    threads: 2
    resources:
             runtime_hrs= 1,
             mem_total_mb = lambda wildcards, attempt: 1024 + 1024 * attempt,
             mem_per_cpu_mb = lambda wildcards, attempt: (1024 + 1024 * attempt) // 2,
    conda:
         '../environment/conda/conda_pyscript.yml'
    params:
          script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
         '{params.script_exec} --debug --input-files {input.fasta} '
         '--output {output.dump} --summary-output {output.summary} '
         '--num-cpu {threads} --genome-size-file {input.faidx} &> {log}'


rule compute_statistics_complete_input_bam:
    input:
        bam = 'input/bam/{sample}.pbn.bam',
        faidx = 'references/assemblies/' + config['use_genome_size'] +'.fasta.fai',
    output:
        dump = 'output/statistics/stat_dumps/{sample}.pbn.bam.pck',
        summary = 'input/bam/{sample}.stats',
    log: 'log/output/statistics/stat_dumps/{sample}.pbn.bam.log',
    benchmark: 'run/output/statistics/stat_dumps/{sample}.pbn.bam.t2.rsrc'
    threads: 2
    resources:
        runtime_hrs= 23,
        mem_total_mb = 4096,
        mem_per_cpu_mb = 2048,
    conda:
         '../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
        '{params.script_exec} --debug --input-files {input.bam} '
        '--output {output.dump} --summary-output {output.summary} '
        '--copy-stats-dump '
        'output/statistics/stat_dumps/{wildcards.sample}.fasta.pck '
        'output/statistics/stat_dumps/{wildcards.sample}.fastq.pck '
        '--copy-summary '
        'input/fasta/{wildcards.sample}.stats '
        'input/fastq/{wildcards.sample}.stats '
        '--num-cpu {threads} --genome-size-file {input.faidx} &> {log}'


rule plot_input_long_reads_statistics:
    input:
        'output/statistics/stat_dumps/{sample}.{file_ext}.pck'
    output:
        'output/plotting/statistics/input_reads/{sample}.{file_ext}.stats.pdf'
    conda:
         '../environment/conda/conda_pyscript.yml'
    priority: 200
    params:
        script_exec = lambda wildcards: find_script_path('plot_sample_stats.py'),
        lower_bound = 6000,
        upper_bound = lambda wildcards: {'ccs': 25000, 'clr': 100000, 'ul': 100000}[wildcards.sample.split('_')[2].split('-')[-1]],
        step_size = lambda wildcards: {'ccs': 500, 'clr': 1000, 'ul': 1000}[wildcards.sample.split('_')[2].split('-')[-1]]
    shell:
        '{params.script_exec} '
        '--pck-input {input} --text-size 11 '
        '--sample-name {wildcards.sample} '
        '--lowest-bin {params.lower_bound} '
        '--highest-bin {params.upper_bound} '
        '--step-size {params.step_size} '
        '--output {output} '


def collect_tag_lists(wildcards, glob_collect=False):
    """
    :param wildcards:
    :return:
    """
    import os

    source_path = os.path.join('output',
                               PATH_STRANDSEQ_DGA_SPLIT,
                               'draft/haplotags/{hap_reads}.{sequence}.tags.{tag_type}.tsv')

    if glob_collect:
        import glob
        pattern = source_path.replace('{sequence}', '*')
        pattern = pattern.format(**dict(wildcards))
        seq_files = glob.glob(pattern)

        if not seq_files:
            raise RuntimeError('collect_tag_lists: no files collected with pattern {}'.format(pattern))

    else:
        reference_folder = os.path.join('output/reference_assembly/clustered', wildcards.sts_reads)
        seq_output_dir = checkpoints.create_assembly_sequence_files.get(folder_path=reference_folder,
                                                                        reference=wildcards.reference).output[0]
        checkpoint_wildcards = glob_wildcards(os.path.join(seq_output_dir, '{sequence}.seq'))

        seq_files = expand(source_path,
                           var_caller=wildcards.var_caller,
                           gq=wildcards.gq,
                           qual=wildcards.qual,
                           reference=wildcards.reference,
                           vc_reads=wildcards.vc_reads,
                           sts_reads=wildcards.sts_reads,
                           hap_reads=wildcards.hap_reads,
                           sequence=checkpoint_wildcards.sequence,
                           tag_type=wildcards.tag_type)
    return seq_files


rule summarize_tagging_splitting_statistics:
    input:
        tags = collect_tag_lists
    output:
        'output/statistics/tag_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.tags.{tag_type}.tsv'
    benchmark:
        'run/output/statistics/tag_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.tags.{tag_type}.rsrc'
    priority: 200
    run:
        try:
            validate_checkpoint_output(input.tags)
            tag_files = input.tags
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{}\n'.format(str(error)))
            tag_files = collect_tag_lists(wildcards, glob_collect=True)

        import os
        import collections as col

        hapcount = col.Counter()
        line_num = 0

        for tsv in sorted(tag_files):
            with open(tsv, 'r') as table:
                for line in table:
                    if line.startswith('#'):
                        continue
                    line_num += 1
                    try:
                        hap = line.split()[1]
                    except IndexError:
                        raise IndexError('Malformed line {} from file {}: "{}"'.format(line_num, os.path.basename(tsv), line))
                    hapcount[hap] += 1

        entry_types = tuple(sorted(hapcount.keys()))
        if entry_types != ('H1', 'H2', 'none'):
            raise ValueError('Haplotype entries in tag list not as expected (H1/H2/none): {}'.format(entry_types))

        total_reads = sum(list(hapcount.values()))

        percent_h1 = round((hapcount['H1'] / total_reads) * 100, 2)
        percent_h2 = round((hapcount['H2'] / total_reads) * 100, 2)
        percent_tag = round(percent_h1 + percent_h2, 2)
        percent_un = round((hapcount['none'] / total_reads) * 100, 2)

        with open(output[0], 'w') as summary:
            _ = summary.write('num_reads\t{}\n'.format(total_reads))
            _ = summary.write('num_hap1\t{}\n'.format(hapcount['H1']))
            _ = summary.write('num_hap2\t{}\n'.format(hapcount['H2']))
            _ = summary.write('num_untagged\t{}\n'.format(hapcount['none']))
            _ = summary.write('percent_hap1\t{}\n'.format(percent_h1))
            _ = summary.write('percent_hap2\t{}\n'.format(percent_h2))
            _ = summary.write('percent_tagged\t{}\n'.format(percent_tag))
            _ = summary.write('percent_untagged\t{}\n'.format(percent_un))
        
