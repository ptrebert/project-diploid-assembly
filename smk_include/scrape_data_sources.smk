
localrules: master_scrape_data_sources

rule master_scrape_data_sources:
    input:
        'input/data_sources/hgsvc_HG00514_pacbio.json',
        'input/data_sources/hgsvc_HG00512_pacbio.json',
        'input/data_sources/hgsvc_HG00513_pacbio.json',
        'input/data_sources/hgsvc_pur-trio_pacbio.json',
        'input/data_sources/hgsvc_NA19240_pacbio.json',
        'input/data_sources/hgsvc_NA19239_pacbio.json',
        'input/data_sources/hgsvc_NA19238_pacbio.json',
        'input/data_sources/hgsvc_JAX_pacbio.json',


rule collect_remote_hgsvc_HG00514_pacbio:
    output:
        'input/data_sources/hgsvc_HG00514_pacbio.json'
    params:
        script_dir = config['script_dir'],
        server = 'ftp.1000genomes.ebi.ac.uk',
        remote_path = 'vol1/ftp/data_collections/HGSVC2/working/20190508_HG00514_PacBioSequel2/',
        collect = ' fastq.gz bam ',
        sort = ' input/fastq/partial/parts input/bam/partial/parts ',
        bam_format = ' --assume-pacbio-native ',
        clr_subreads = ' --assume-clr-subreads',
        file_infix = ' hgsvc_pbsq2- '
    log:
        'log/input/data_sources/hgsvc_chs_HG00514_pacbio.log'
    conda:
        '../environment/conda/conda_pyscript.yml'
    shell:
        '{params.script_dir}/scan_remote_path.py --debug '
            ' --server {params.server} --ftp-path {params.remote_path} '
            ' --collect-files {params.collect} --sort-files {params.sort} '
            ' {params.bam_format} {params.clr_subreads} --file-infix {params.file_infix}'
            ' --output {output} &> {log}'


rule collect_remote_hgsvc_HG00512_pacbio:
    output:
        'input/data_sources/hgsvc_HG00512_pacbio.json'
    params:
        script_dir = config['script_dir'],
        server = 'ftp.1000genomes.ebi.ac.uk',
        remote_path = 'vol1/ftp/data_collections/HGSVC2/working/20191031_CHS_PacBio_HG00512_HiFi/',
        collect = ' fastq.gz bam ',
        sort = ' input/fastq/partial/parts input/bam/partial/parts ',
        bam_format = ' --assume-pacbio-native ',
        file_infix = ' hgsvc_pbsq2- '
    log:
        'log/input/data_sources/hgsvc_chs_HG00512_pacbio.log'
    conda:
        '../environment/conda/conda_pyscript.yml'
    shell:
        '{params.script_dir}/scan_remote_path.py --debug '
            ' --server {params.server} --ftp-path {params.remote_path} '
            ' --collect-files {params.collect} --sort-files {params.sort} '
            ' {params.bam_format} --file-infix {params.file_infix}'
            ' --output {output} &> {log}'


rule collect_remote_hgsvc_HG00513_pacbio:
    output:
        'input/data_sources/hgsvc_HG00513_pacbio.json'
    params:
        script_dir = config['script_dir'],
        server = 'ftp.1000genomes.ebi.ac.uk',
        remote_path = 'vol1/ftp/data_collections/HGSVC2/working/20191031_CHS_PacBio_HG00513_HiFi/',
        collect = ' fastq.gz bam ',
        sort = ' input/fastq/partial/parts input/bam/partial/parts ',
        bam_format = ' --assume-pacbio-native ',
        file_infix = ' hgsvc_pbsq2- '
    log:
        'log/input/data_sources/hgsvc_chs_HG00513_pacbio.log'
    conda:
        '../environment/conda/conda_pyscript.yml'
    shell:
        '{params.script_dir}/scan_remote_path.py --debug '
            ' --server {params.server} --ftp-path {params.remote_path} '
            ' --collect-files {params.collect} --sort-files {params.sort} '
            ' {params.bam_format} --file-infix {params.file_infix}'
            ' --output {output} &> {log}'


rule collect_remote_hgsvc_pur_trio_pacbio:
    output:
        'input/data_sources/hgsvc_pur-trio_pacbio.json'
    params:
        script_dir = config['script_dir'],
        server = 'ftp.1000genomes.ebi.ac.uk',
        remote_path = 'vol1/ftp/data_collections/HGSVC2/working/20190925_PUR_PacBio_HiFi/',
        collect = ' fastq.gz bam ',
        sort = ' input/fastq/partial/parts input/bam/partial/parts ',
        bam_format = ' --assume-pacbio-native ',
        file_infix = ' hgsvc_pbsq2- '
    log:
        'log/input/data_sources/hgsvc_pur-trio_pacbio.log'
    conda:
        '../environment/conda/conda_pyscript.yml'
    shell:
        '{params.script_dir}/scan_remote_path.py --debug '
            ' --server {params.server} --ftp-path {params.remote_path} '
            ' --collect-files {params.collect} --sort-files {params.sort} '
            ' {params.bam_format} --file-infix {params.file_infix}'
            ' --output {output} &> {log}'


rule collect_remote_hgsvc_NA19240_pacbio:
    output:
        'input/data_sources/hgsvc_NA19240_pacbio.json'
    params:
        script_dir = config['script_dir'],
        server = 'ftp.1000genomes.ebi.ac.uk',
        remote_path = 'vol1/ftp/data_collections/HGSVC2/working/20191005_YRI_PacBio_NA19240_HiFi/',
        collect = ' fastq.gz bam ',
        sort = ' input/fastq/partial/parts input/bam/partial/parts ',
        bam_format = ' --assume-pacbio-native ',
        file_infix = ' hgsvc_pbsq2- '
    log:
        'log/input/data_sources/hgsvc_yri_NA19240_pacbio.log'
    conda:
        '../environment/conda/conda_pyscript.yml'
    shell:
        '{params.script_dir}/scan_remote_path.py --debug '
            ' --server {params.server} --ftp-path {params.remote_path} '
            ' --collect-files {params.collect} --sort-files {params.sort} '
            ' {params.bam_format} --file-infix {params.file_infix}'
            ' --output {output} &> {log}'


rule collect_remote_hgsvc_NA19238_pacbio:
    output:
        'input/data_sources/hgsvc_NA19238_pacbio.json'
    params:
        script_dir = config['script_dir'],
        server = 'ftp.1000genomes.ebi.ac.uk',
        remote_path = 'vol1/ftp/data_collections/HGSVC2/working/20191205_YRI_PacBio_NA19238_HIFI/',
        collect = ' fastq.gz bam ',
        sort = ' input/fastq/partial/parts input/bam/partial/parts ',
        bam_format = ' --assume-pacbio-native ',
        file_infix = ' hgsvc_pbsq2- '
    log:
        'log/input/data_sources/hgsvc_yri_NA19238_pacbio.log'
    conda:
        '../environment/conda/conda_pyscript.yml'
    shell:
        '{params.script_dir}/scan_remote_path.py --debug '
            ' --server {params.server} --ftp-path {params.remote_path} '
            ' --collect-files {params.collect} --sort-files {params.sort} '
            ' {params.bam_format} --file-infix {params.file_infix}'
            ' --output {output} &> {log}'


rule collect_remote_hgsvc_NA19239_pacbio:
    output:
        'input/data_sources/hgsvc_NA19239_pacbio.json'
    params:
        script_dir = config['script_dir'],
        server = 'ftp.1000genomes.ebi.ac.uk',
        remote_path = 'vol1/ftp/data_collections/HGSVC2/working/20191205_YRI_PacBio_NA19239_HIFI/',
        collect = ' fastq.gz bam ',
        sort = ' input/fastq/partial/parts input/bam/partial/parts ',
        bam_format = ' --assume-pacbio-native ',
        file_infix = ' hgsvc_pbsq2- '
    log:
        'log/input/data_sources/hgsvc_yri_NA19239_pacbio.log'
    conda:
        '../environment/conda/conda_pyscript.yml'
    shell:
        '{params.script_dir}/scan_remote_path.py --debug '
            ' --server {params.server} --ftp-path {params.remote_path} '
            ' --collect-files {params.collect} --sort-files {params.sort} '
            ' {params.bam_format} --file-infix {params.file_infix}'
            ' --output {output} &> {log}'


rule collect_remote_hgsvc_jax_pacbio:
    output:
        'input/data_sources/hgsvc_JAX_pacbio.json'
    params:
        script_dir = config['script_dir'],
        server = 'ftp.1000genomes.ebi.ac.uk',
        remote_path = 'vol1/ftp/data_collections/HGSVC2/working/20200108_PacBio_CLR_JAX/',
        collect = ' bam ',
        sort = ' input/bam/partial/parts ',
        bam_format = ' --assume-pacbio-native ',
        clr_subreads = ' --assume-clr-subreads',
        file_infix = ' hgsvc_pbsq2- '
    log:
        'log/input/data_sources/hgsvc_JAX_pacbio.log'
    conda:
        '../environment/conda/conda_pyscript.yml'
    shell:
        '{params.script_dir}/scan_remote_path.py --debug '
            ' --server {params.server} --ftp-path {params.remote_path} '
            ' --collect-files {params.collect} --sort-files {params.sort} '
            ' {params.bam_format} {params.clr_subreads} --file-infix {params.file_infix}'
            ' --output {output} &> {log}'