
localrules: create_cluster_log_folders


rule create_cluster_log_folders:
    output:
        touch('output/check_files/environment/cluster_log_folders.ok')
    shell:
        'mkdir -p log/cluster_jobs/stderr && mkdir -p log/cluster_jobs/stdout'


rule create_conda_environment_shell_tools:
    output:
        'output/check_files/environment/conda_shelltools.ok'
    log:
        'log/output/check_files/environment/conda_shelltools.log'
    params:
        script_exec = lambda wildcards: find_script_path('inspect_environment.py', 'utilities')
    conda:
        '../environment/conda/conda_shelltools.yml'
    shell:
         '{params.script_exec} '
         '--export-conda-env --outfile {output} --logfile {log}'


rule create_conda_environment_pacbio_tools:
    output:
        'output/check_files/environment/conda_pbtools.ok'
    log:
        'log/output/check_files/environment/conda_pbtools.log'
    params:
        script_exec = lambda wildcards: find_script_path('inspect_environment.py', 'utilities')
    conda:
        '../environment/conda/conda_pbtools.yml'
    shell:
         '{params.script_exec} '
         '--export-conda-env --outfile {output} --logfile {log}'


rule create_conda_environment_r_script:
    """
    Conda environment for execution of R scripts
    SaaRclust / breakpointR / StrandPhaseR
    """
    output:
        'output/check_files/environment/conda_rscript.ok'
    log:
        'log/output/check_files/environment/conda_rscript.log'
    params:
        script_exec = lambda wildcards: find_script_path('inspect_environment.py', 'utilities')
    conda:
        '../environment/conda/conda_rscript.yml'
    shell:
         '{params.script_exec} '
         '--export-conda-env --outfile {output} --logfile {log}'


rule create_conda_environment_r_tools:
    """
    Conda environment for execution of tools with
    an R dependency such as Quast/BUSCO
    """
    output:
        'output/check_files/environment/conda_rtools.ok'
    log:
        'log/output/check_files/environment/conda_rtools.log'
    params:
        script_exec = lambda wildcards: find_script_path('inspect_environment.py', 'utilities')
    conda:
        '../environment/conda/conda_rtools.yml'
    shell:
         '{params.script_exec} '
         '--export-conda-env --outfile {output} --logfile {log}'


rule create_conda_environment_bio_tools:
    output:
        'output/check_files/environment/conda_biotools.ok'
    log:
        'log/output/check_files/environment/conda_biotools.log'
    params:
        script_exec = lambda wildcards: find_script_path('inspect_environment.py', 'utilities')
    conda:
        '../environment/conda/conda_biotools.yml'
    shell:
         '{params.script_exec} '
         '--export-conda-env --outfile {output} --logfile {log}'


rule create_conda_environment_pyscript:
    output:
        'output/check_files/environment/conda_pyscript.ok'
    log:
        'log/output/check_files/environment/conda_pyscript.log'
    params:
        script_exec = lambda wildcards: find_script_path('inspect_environment.py', 'utilities')
    conda:
        '../environment/conda/conda_pyscript.yml'
    shell:
         '{params.script_exec} '
         '--export-conda-env --outfile {output} --logfile {log}'


rule inspect_hpc_module_singularity:
    output:
        'output/check_files/environment/module_singularity.ok'
    log:
        'log/output/check_files/environment/module_singularity.log'
    params:
        script_exec = lambda wildcards: find_script_path('inspect_environment.py', 'utilities'),
        singularity_module = config['env_module_singularity']
    #envmodules:
    #    config['env_module_singularity']
    shell:
         '{params.script_exec} '
         '--outfile {output} --logfile {log}'


rule check_singularity_version:
    input:
        'output/check_files/environment/module_singularity.ok'
    output:
        'output/check_files/environment/singularity_version.ok'
    log:
        'log/output/check_files/environment/singularity_version.log'
    #envmodules:
    #    config['env_module_singularity']
    params:
        script_exec = lambda wildcards: find_script_path('version_checker.py', 'utilities'),
        min_version = '3.1.0',  # due to container format change between v2 and v3
    shell:
        'singularity --version | '
        '{params.script_exec} '
        '--outfile {output} --logfile {log} '
        '--at-least {params.min_version}'


rule download_shasta_executable:
    input:
        'output/check_files/environment/conda_biotools.ok'
    output:
        'output/check_files/environment/shasta_version.ok'
    log:
        'log/output/check_files/environment/shasta_version.log'
    conda:
        '../environment/conda/conda_biotools.yml'
    params:
        shasta_url = "https://github.com/chanzuckerberg/shasta/releases/download/{version}/shasta-Linux-{version}".format(**{'version': config['shasta_version']}),
        shasta_ver = config['shasta_version'],
        script_exec = lambda wildcards: find_script_path('downloader.py', 'utilities'),
    threads: 2
    resources:
        mem_total_mb = 2048,
        mem_per_cpu_mb = 1024
    shell:
         '{params.script_exec} --debug '
            '--shasta-exec {params.shasta_url} '
            '--shasta-version {params.shasta_ver} '
            '--shasta-path $CONDA_PREFIX/bin/shasta '
            '--output {output} '
            '&> {log}'


rule install_rlib_saarclust:
    input:
        rules.create_conda_environment_r_script.output
    output:
         check = touch('output/check_files/R_setup/saarclust_ver-{}.ok'.format(config['git_commit_saarclust']))
    log:
        'log/output/check_files/R_setup/saarclust_ver-{}.log'.format(config['git_commit_saarclust'])
    conda:
        '../environment/conda/conda_rscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('install_saarclust.R'),
        version = config['git_commit_saarclust']
    resources:
        mem_total_mb = 4096,
        mem_per_cpu_mb = 4096
    shell:
        'TAR=$(which tar) {params.script_exec} {params.version} &> {log}'


if int(config['git_commit_version']) < 8:
    INSTALL_BREAKPOINTR_OUTPUT = 'output/check_files/R_setup/breakpointr.ok'
    BREAKPOINTR_VERSION = -1
else:
    INSTALL_BREAKPOINTR_OUTPUT = 'output/check_files/R_setup/breakpointr_ver-{}.ok'.format(config['git_commit_breakpointr'])
    BREAKPOINTR_VERSION = config['git_commit_breakpointr']


rule install_rlib_breakpointr:
    input:
        rules.install_rlib_saarclust.output.check
    output:
         check = touch(INSTALL_BREAKPOINTR_OUTPUT)
    log:
        'log/output/check_files/R_setup/breakpointr.log'
    conda:
        '../environment/conda/conda_rscript.yml'
    resources:
        mem_total_mb = 4096,
        mem_per_cpu_mb = 4096
    params:
        script_exec = lambda wildcards: find_script_path('install_breakpointr.R'),
        version = BREAKPOINTR_VERSION
    shell:
        'TAR=$(which tar) {params.script_exec} {params.version} &> {log}'


rule install_rlib_strandphaser:
    input:
        rules.install_rlib_breakpointr.output.check
    output:
         check = touch('output/check_files/R_setup/strandphaser_ver-{}.ok'.format(config['git_commit_strandphaser']))
    log:
        'log/output/check_files/R_setup/strandphaser_ver-{}.log'.format(config['git_commit_strandphaser'])
    conda:
        '../environment/conda/conda_rscript.yml'
    resources:
        mem_total_mb = 4096,
        mem_per_cpu_mb = 4096
    params:
        script_exec = lambda wildcards: find_script_path('install_strandphaser.R'),
        version = config['git_commit_strandphaser']
    shell:
        'TAR=$(which tar) {params.script_exec} {params.version} &> {log}'


rule download_quast_busco_databases:
    input:
        rules.create_conda_environment_r_tools.output
    output:
        'output/check_files/quast-lg/busco_db_download.ok'
    conda:
        '../environment/conda/conda_rtools.yml'
    shell:
        'quast-download-busco &> {output}'