

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


rule check_singularity_version:
    output:
        'output/check_files/environment/singularity_version.ok'
    envmodules:
        config['env_module_singularity']
    run:
        import subprocess as sp

        try:
            sing_ver = sp.check_output('singularity --version',
                                        stderr=sp.STDOUT,
                                        shell=True)
            sing_ver = sing_ver.decode('utf-8')
            version_string = sing_ver.strip().split()[-1]
            major, minor = version_string.split('.')[:2]
            if int(major) >= 3 and int(minor) > 0:
                with open(output[0], 'w') as dump:
                    _ = dump.write(sing_ver)
            else:
                raise ValueError('Incompatible Singularity version (>3.0 required): {}'.format(sing_ver))
        except sp.CalledProcessError as spe:
            rc = spe.returncode
            err_msg = str(spe.output)
            raise ValueError('Could not determine Singularity version (>3.0 required): {} / {}'.format(rc, err_msg))


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