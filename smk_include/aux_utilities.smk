
rule compute_md5_checksum:
    input:
        '{filepath}'
    output:
        '{filepath}.md5'
    shell:
        "md5sum {input} > {output}"


rule samtools_index_cram_alignment:
    input:
        cram = '{filepath}.cram'
    output:
        crai = '{filepath}.cram.crai'
    threads: config['num_cpu_low']
    shell:
        "samtools index -@ {threads} {input.cram}"


rule samtools_index_bam_alignment:
    input:
        bam = '{filepath}'
    output:
        bai = '{filepath}.bai'
    threads: config['num_cpu_low']
    shell:
        "samtools index -@ {threads} {input.bam}"


rule pb_bamtools_index_bam_alignment:
    input:
        pbn_bam = '{filepath}.pbn.bam'
    output:
        pbi = '{filepath}.pbn.bam.pbi'
    benchmark:
        'run/{filepath}.create-pbi.rsrc'
    resources:
        runtime_hrs = lambda wildcards: 0 if '.part' in wildcards.filepath else 1
    conda:
        config['conda_env_pbtools']
    shell:
        'pbindex {input}'


rule pb_bam2x_dump_fastq:
    input:
        pbn_bam = 'input/bam/{folder_path}/{pbn_sample}{sample_type}.pbn.bam',
        pbn_idx = 'input/bam/{folder_path}/{pbn_sample}{sample_type}.pbn.bam.pbi'
    output:
        'input/fastq/{folder_path}/{pbn_sample}{sample_type}.fastq.gz'
    log:
        'log/input/bam/{folder_path}/{pbn_sample}{sample_type}.dump.log'
    benchmark:
        'run/input/bam/{folder_path}/{pbn_sample}{sample_type}.dump.rsrc'
    wildcard_constraints:
        pbn_sample = '(' + '|'.join(config['partial_pbn_samples'] + config['complete_pbn_samples']) + ')',
        sample_type = '(_[0-9x]+|\.part[0-9]+)'
    resources:
        runtime_hrs = lambda wildcards, input: 1 if '.part' in input.pbn_bam else 6
    conda:
        config['conda_env_pbtools']
    params:
        out_prefix = lambda wildcards, output: output[0].rsplit('.', 2)[0]
    shell:
        'bam2fastq -c 5 -o {params.out_prefix} {input.pbn_bam} &> {log}'



rule samtools_convert_sam_to_bam:
    input:
        sam = '{folder_path}.sam'
    output:
        bam = '{folder_path}.sam.bam'
    threads: config['num_cpu_low']
    benchmark:
        'run/{folder_path}.sam-convert.rsrc'
    shell:
        "samtools view -o {output.bam} -b -@ {threads} {input.sam}"


rule samtools_position_sort_bam_alignment:
    input:
        unsorted_bam = '{folder_path}.sam.bam'
    output:
        sorted_bam = '{folder_path}.psort.sam.bam'
    threads: config['num_cpu_low']
    benchmark:
        'run/{folder_path}.psort-bam.rsrc'
    resources:
        mem_per_cpu_mb = 20 * 1024,  # 20G per CPU
        mem_total_mb = config['num_cpu_low'] * (20 * 1024)
    params:
        mem_per_thread = '20G'
    shell:
        'samtools sort -m {params.mem_per_thread} --threads {threads} ' \
            '-o {output.sorted_bam} {input.unsorted_bam}'


rule samtools_index_fasta:
    input:
        '{filepath}.fasta'
    output:
        '{filepath}.fasta.fai'
    shell:
        'samtools faidx {input}'


rule bgzip_file_copy:
    input:
        '{filepath}'
    output:
        '{filepath}.bgz'
    shell:
        'bgzip -c {input} > {output}'


rule bcftools_index_bgzipped_file:
    input:
        '{filepath}.bgz'
    output:
        '{filepath}.bgz.tbi'
    shell:
        'bcftools index --tbi {input}'


checkpoint create_assembly_sequence_files:
    input:
        '{folder_path}/{reference}.fasta.fai'
    output:
        directory('{folder_path}/{reference}/sequences')
    run:
        output_dir = output[0]
        os.makedirs(output_dir, exist_ok=True)
        with open(input[0], 'r') as fai:
            for line in fai:
                seq_name = line.split('\t')[0]
                output_path = os.path.join(output_dir, seq_name + '.seq')
                with open(output_path, 'w') as dump:
                    _ = dump.write(line)


rule generate_bwa_index:
    input:
        reference = '{folder_path}/{reference}.fasta'
    output:
        '{folder_path}/{reference}/bwa_index/{reference}.amb',
        '{folder_path}/{reference}/bwa_index/{reference}.ann',
        '{folder_path}/{reference}/bwa_index/{reference}.bwt',
        '{folder_path}/{reference}/bwa_index/{reference}.pac',
        '{folder_path}/{reference}/bwa_index/{reference}.sa'
    log:
        'log/{folder_path}/bwa_index/{reference}.log'
    benchmark:
        'run/{folder_path}/bwa_index/{reference}.rsrc'
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 6144
    params:
        prefix = lambda wildcards, output: output[0].split('.')[0]
    shell:
        'bwa index -p {params.prefix} {input.reference} &> {log}'


rule no_singularity_mock_output:
    output:
        touch('output/check_files/environment/singularity.false')


rule check_singularity_version:
    output:
        'output/check_files/environment/singularity_version.ok'
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


rule singularity_pull_container:
    input:
        'output/check_files/environment/singularity_version.ok'
    output:
        'output/container/{hub}/{repo}/{tool}_{version}.sif'
    log:
        'log/output/container/{hub}/{repo}/{tool}_{version}.pull.log'
    params:
        pull_folder = lambda wildcards: os.path.join(os.getcwd(), 'output', 'container', wildcards.hub, wildcards.repo)
    shell:
        'SINGULARITY_PULLFOLDER={params.pull_folder} singularity pull ' \
            '{wildcards.hub}://{wildcards.repo}/{wildcards.tool}:{wildcards.version} &> {log}'


def collect_strandseq_alignments(wildcards):
    """
    """
    individual, project, platform_spec = wildcards.sts_reads.split('_')[:3]
    platform, spec = platform_spec.split('-')

    requests_dir = checkpoints.create_bioproject_download_requests.get(sts_reads=wildcards.sts_reads).output[0]

    search_path = os.path.join(requests_dir, '{individual}_{project}_{platform_spec}_{lib_id}_{run_id}_1.request')

    checkpoint_wildcards = glob_wildcards(search_path)

    bam_files = expand(
        'output/alignments/strandseq_to_reference/{reference}/{sts_reads}/{individual}_{project}_{platform}-npe_{lib_id}.mrg.psort.mdup.sam.bam{ext}',
        reference=wildcards.reference,
        individual=individual,
        sts_reads=wildcards.sts_reads,
        project=project,
        platform=platform,
        lib_id=checkpoint_wildcards.lib_id,
        ext=['', '.bai'])
    return sorted(bam_files)


def load_preset_file(wildcards, input):
    """
    Save load for parameters from files that do not
    exist during a dry run, leading to a FileNotFoundError
    raised by Snakemake. Replaces construct

    params:
        preset = lambda wildcards, input: open(input.preset).read().strip()

    with

    params:
        preset = load_preset_file

    """
    if not hasattr(input, 'preset'):
        raise AttributeError('Input does not have a "preset" attribute: {}'.format(input))
    file_path = input.preset
    if not os.path.isfile(file_path):
        preset = 'PRESET-DRY-RUN'
    else:
        with open(file_path, 'r') as dump:
            preset = dump.read().strip()
            assert preset, 'Empty preset file: {}'.format(file_path)
    return preset


def load_fofn_file(input, prefix='', sep=' '):
    """
    Save load for list of filenames from a file that
    does not exist during a dry run, leading to a
    FileNotFoundError raised by Snakemake.

    Replaces construct

    params:
        preset = lambda wildcards, input: open(input.fofn).read().strip()

    with

    params:
        preset = lambda wildcards, input: load_fofn_file(input, prefix_string, separator_string)
    """
    if not hasattr(input, 'fofn'):
        raise AttributeError('Input does not have "fofn" attribute: {}'.format(input))
    file_path = input.fofn
    if not os.path.isfile(file_path):
        file_list = ['FOFN-DRY-RUN']
    else:
        with open(file_path, 'r') as dump:
            file_list = sorted([l.strip() for l in dump.readlines()])
            assert file_list, 'Empty fofn file: {}'.format(file_path)
    file_list = prefix + sep.join(file_list)
    return file_list


def load_seq_length_file(wildcards, input):
    """
    Follows same idea as two function above
    """
    if not hasattr(input, 'seq_info'):
        raise AttributeError('Input does not have "seq_info" attribute: {}'.format(input))
    file_path = input.seq_info
    if not os.path.isfile(file_path):
        seq_len = 'SEQLEN-DRY-RUN'
    else:
        with open(file_path, 'r') as info:
            seq_len = info.readline().split()[1]
            try:
                _ = int(seq_len)
            except ValueError:
                raise AssertionError('Extracted seq. length is not an integer: {} / {}'.format(seq_len, file_path))
    return seq_len