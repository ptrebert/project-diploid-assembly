
include: '../aux_utilities.smk'

rule create_conda_environment_compile:
    output:
        'output/check_files/environment/conda_compile.ok'
    log:
        'log/output/check_files/environment/conda_compile.log'
    params:
        script_exec = lambda wildcards: find_script_path('inspect_environment.py', 'utilities')
    conda:
        '../../environment/conda/conda_compile.yml'
    shell:
        '{params.script_exec} '
        '--export-conda-env --outfile {output} --logfile {log}'


rule install_source_bifrost:
    """
    Bioconda version of Bifrost is built with AVX2 enabled,
    which does not work on older machines. Additionally, recent
    bug fix for graph querying has not been part of any release:
    github.com/pmelsted/bifrost/issues/20
    github.com/pmelsted/bifrost/issues/21
    """
    input:
        'output/check_files/environment/conda_compile.ok'
    output:
        touch('output/check_files/src_build/install_bifrost.ok')
    log:
       'output/check_files/src_build/install_bifrost.ok'
    conda:
        '../../environment/conda/conda_compile.yml'
    params:
        repo_folder = 'output/repositories'
    shell:
         'rm -rf {params.repo_folder}/bifrost && '
         'mkdir -p {params.repo_folder} && '
         'cd {params.repo_folder} && '
         'git clone https://github.com/pmelsted/bifrost.git && '
         'cd bifrost && '
         'git checkout f1e48443f3576429590b65d24b998c88a52fb6d4 && '
         'mkdir build && cd build && '
         'cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX .. && '
         'make && make install &> {log}'