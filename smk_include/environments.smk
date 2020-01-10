

rule create_conda_environment_shell_tools:
    output:
        'output/check_files/environment/conda_shelltools.ok'
    log:
        'log/output/check_files/environment/conda_shelltools.log'
    params:
        script_dir = config['script_dir']
    conda:
        '../environment/conda/conda_shelltools.yml'
    shell:
         '{params.script_dir}/utilities/inspect_environment.py --outfile {output} --logfile {log}'


rule create_conda_environment_pacbio_tools:
    output:
        'output/check_files/environment/conda_pbtools.ok'
    log:
        'log/output/check_files/environment/conda_pbtools.log'
    params:
        script_dir = config['script_dir']
    conda:
        '../environment/conda/conda_pbtools.yml'
    shell:
         '{params.script_dir}/utilities/inspect_environment.py --outfile {output} --logfile {log}'


rule create_conda_environment_r_tools:
    output:
        'output/check_files/environment/conda_rtools.ok'
    log:
        'log/output/check_files/environment/conda_rtools.log'
    params:
        script_dir = config['script_dir']
    conda:
        '../environment/conda/conda_rtools.yml'
    shell:
         '{params.script_dir}/utilities/inspect_environment.py --outfile {output} --logfile {log}'


rule create_conda_environment_bio_tools:
    output:
        'output/check_files/environment/conda_biotools.ok'
    log:
        'log/output/check_files/environment/conda_biotools.log'
    params:
        script_dir = config['script_dir']
    conda:
        '../environment/conda/conda_biotools.yml'
    shell:
         '{params.script_dir}/utilities/inspect_environment.py --outfile {output} --logfile {log}'


rule create_conda_environment_pyscript:
    output:
        'output/check_files/environment/conda_pyscript.ok'
    log:
        'log/output/check_files/environment/conda_pyscript.log'
    params:
        script_dir = config['script_dir']
    conda:
        '../environment/conda/conda_pyscript.yml'
    shell:
         '{params.script_dir}/utilities/inspect_environment.py --outfile {output} --logfile {log}'


rule inspect_hpc_module_singularity:
    output:
        'output/check_files/environment/module_singularity.ok'
    log:
        'log/output/check_files/environment/module_singularity.ok'
    params:
        script_dir = config['script_dir'],
        singularity_module = config['env_module_singularity']
    #envmodules:
    #    config['env_module_singularity']
    shell:
#         'module load {params.singularity_module} ; '
         '{params.script_dir}/utilities/inspect_environment.py --outfile {output} --logfile {log} ; '
#         'module unload {params.singularity_module} '


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
        script_dir = config['script_dir'],
        min_version = '3.1.0',  # due to container format change between v2 and v3
        singularity_module = config['env_module_singularity']
    shell:
#        'module load {params.singularity_module} ; '
        'singularity --version | '
        '{params.script_dir}/utilities/version_checker.py '
        '--outfile {output} --logfile {log} '
        '--at-least {params.min_version} ; '
#        'module unload {params.singularity_module}'


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
        script_dir = config['script_dir']
    threads: 2
    resources:
        mem_total_mb = 2048,
        mem_per_cpu_mb = 1024
    shell:
         '{params.script_dir}/utilities/downloader.py --debug '
            '--shasta-exec {params.shasta_url} '
            '--shasta-version {params.shasta_ver} '
            '--shasta-path $CONDA_PREFIX/bin/shasta '
            '--output {output} '
            '&> {log}'


rule install_rlib_saarclust:
    input:
        'output/check_files/environment/conda_rtools.ok'
    output:
         check = touch('output/check_files/R_setup/saarclust_ver-{}.ok'.format(config['git_commit_saarclust']))
    log:
        'log/output/check_files/R_setup/saarclust_ver-{}.log'.format(config['git_commit_saarclust'])
    conda:
        '../environment/conda/conda_rtools.yml'
    params:
        script_dir = config['script_dir'],
        version = config['git_commit_saarclust']
    resources:
        mem_total_mb = 4096,
        mem_per_cpu_mb = 4096
    shell:
        'TAR=$(which tar) {params.script_dir}/install_saarclust.R {params.version} &> {log}'


rule install_rlib_breakpointr:
    input:
        #'output/check_files/R_setup/saarclust_ver-{}.ok'.format(config['git_commit_saarclust'])
        rules.install_rlib_saarclust.output.check
    output:
         check = touch('output/check_files/R_setup/breakpointr.ok')
    log:
        'log/output/check_files/R_setup/breakpointr.log'
    conda:
        '../environment/conda/conda_rtools.yml'
    resources:
        mem_total_mb = 4096,
        mem_per_cpu_mb = 4096
    params:
        script_dir = config['script_dir']
    shell:
        'TAR=$(which tar) {params.script_dir}/install_breakpointr.R &> {log}'


rule install_rlib_strandphaser:
    input:
        rules.install_rlib_breakpointr.output.check
    output:
         check = touch('output/check_files/R_setup/strandphaser_ver-{}.ok'.format(config['git_commit_strandphaser']))
    log:
        'log/output/check_files/R_setup/strandphaser_ver-{}.log'.format(config['git_commit_strandphaser'])
    conda:
        '../environment/conda/conda_rtools.yml'
    resources:
        mem_total_mb = 4096,
        mem_per_cpu_mb = 4096
    params:
        script_dir = config['script_dir'],
        version = config['git_commit_strandphaser']
    shell:
        'TAR=$(which tar) {params.script_dir}/install_strandphaser.R {params.version} &> {log}'


rule download_quast_busco_databases:
    input:
        rules.install_rlib_strandphaser.output.check
    output:
        'output/check_files/quast-lg/busco_db_download.ok'
    conda:
        '../environment/conda/conda_rtools.yml'
    shell:
        'quast-download-busco &> {output}'