
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


rule plot_saarclust_nhr_assembly_diagnostic_output:
    """
    The default for SaaRclust is to concatenate all individual contigs
    per cluster into a single sequence. During this process, ordering
    information is lost, hence the following output is not part of this rule:
    ordering = 'output/plotting/saarclust_diagnostics/{folder_path}/{reference}_map-to_{aln_reference}.ordering.pdf',
    """
    input:
        setup_ok = 'output/check_files/R_setup/saarclust_ver-{}.ok'.format(config['git_commit_saarclust']),
        ctg_ref_aln = 'output/alignments/contigs_to_reference/reference_assembly/{folder_path}/{reference}_map-to_{aln_reference}.bed'
    output:
        clustering = 'output/plotting/saarclust_diagnostics/reference_assembly/{folder_path}/{reference}_map-to_{aln_reference}.clustering.pdf',
        orienting = 'output/plotting/saarclust_diagnostics/reference_assembly/{folder_path}/{reference}_map-to_{aln_reference}.orienting.pdf',
    log:
       'log/output/plotting/saarclust_diagnostics/reference_assembly/{folder_path}/{reference}_map-to_{aln_reference}.saarclust-diagnostics.log'
    benchmark:
        'run/output/plotting/saarclust_diagnostics/reference_assembly/{folder_path}/{reference}_map-to_{aln_reference}.saarclust-diagnostics.rsrc'
    conda:
        '../environment/conda/conda_rscript.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt
    params:
        script_exec = lambda wildcards: find_script_path('plot_saarclust_diagnostics.R'),
        out_prefix = lambda wildcards: os.path.join(
            'output', 'plotting', 'saarclust_diagnostics', reference_assembly, wildcards.folder_path,
            wildcards.reference + '_map-to_' + wildcards.aln_reference)
    shell:
         '{params.script_exec} {input.ctg_ref_aln} hg38 {params.out_prefix} FALSE &> {log}'


rule plot_saarclust_haploid_assembly_diagnostic_output:
    """
    
    """
    input:
        setup_ok = 'output/check_files/R_setup/saarclust_ver-{}.ok'.format(config['git_commit_saarclust']),
        ctg_ref_aln = 'output/alignments/contigs_to_reference/diploid_assembly/{folder_path}/{file_name}_map-to_{aln_reference}.bed'
    output:
        clustering = 'output/plotting/saarclust_diagnostics/diploid_assembly/{folder_path}/{file_name}_map-to_{aln_reference}.clustering.pdf',
        orienting = 'output/plotting/saarclust_diagnostics/diploid_assembly/{folder_path}/{file_name}_map-to_{aln_reference}.orienting.pdf',
        ordering = 'output/plotting/saarclust_diagnostics/diploid_assembly/{folder_path}/{file_name}_map-to_{aln_reference}.ordering.pdf',
    log:
       'log/output/plotting/saarclust_diagnostics/diploid_assembly/{folder_path}/{file_name}_map-to_{aln_reference}.saarclust-diagnostics.log'
    benchmark:
        'run/output/plotting/saarclust_diagnostics/diploid_assembly/{folder_path}/{file_name}_map-to_{aln_reference}.saarclust-diagnostics.rsrc'
    conda:
        '../environment/conda/conda_rscript.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt
    params:
        script_exec = lambda wildcards: find_script_path('plot_saarclust_diagnostics.R'),
        out_prefix = lambda wildcards: os.path.join(
            'output', 'plotting', 'saarclust_diagnostics', diploid_assembly, wildcards.folder_path,
            wildcards.file_name + '_map-to_' + wildcards.aln_reference)
    shell:
         '{params.script_exec} {input.ctg_ref_aln} hg38 {params.out_prefix} TRUE &> {log}'