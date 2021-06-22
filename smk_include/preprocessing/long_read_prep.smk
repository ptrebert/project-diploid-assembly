import os as os
import fnmatch as fnm

localrules: run_all

HIFI_ROOT_PATH = config.get('hifi_root', 'no_hifi_root_path_set')
if not os.path.isdir(HIFI_ROOT_PATH):
    raise RuntimeError('HiFi root path does not exist: {}'.format(HIFI_ROOT_PATH))


def collect_uncompressed_fastq():

    filepaths = []
    filenames = []

    for root, dirs, files in os.walk(HIFI_ROOT_PATH, followlinks=False):
        fastq_files = fnm.filter(files, '*.fastq')
        filenames.extend([fn.rsplit('.', 1)[0] for fn in fastq_files])
        gz_fastq_files = [os.path.join(root, f + '.gz') for f in fastq_files]
        filepaths.extend(gz_fastq_files)

    constraint = '(' + '|'.join(filenames) + ')'

    return filepaths, constraint


FASTQ_TO_GZIP, CONSTRAINT_FASTQ = collect_uncompressed_fastq()


def collect_bam_only():

    filepaths = []
    filenames = []

    for root, dirs, files in os.walk(HIFI_ROOT_PATH, followlinks=False):
        fastq_files = fnm.filter(files, '*.fastq*')
        if fastq_files:
            # skip over dirs containing FASTQ files
            continue
        bam_files = fnm.filter(files, '*.bam')
        filenames.extend([fn.rsplit('.', 1)[0] for fn in bam_files])
        gz_fastq_files = [os.path.join(root, f.rsplit('.', 1)[0] + '.fastq.gz') for f in fastq_files]
        filepaths.extend(gz_fastq_files)

    constraint = '(' + '|'.join(filenames) + ')'

    return filepaths, constraint


BAM_TO_FASTQ, CONSTRAINT_BAM = collect_bam_only()


rule run_all:
    input:
        FASTQ_TO_GZIP,
        BAM_TO_FASTQ



rule convert_bam_to_fastq:
    input:
        bam = os.path.join(HIFI_ROOT_PATH, '{subfolder}', '{filename}.bam'),
        idx = os.path.join(HIFI_ROOT_PATH, '{subfolder}', '{filename}.bam.pbi')
    output:
        fastq = os.path.join(HIFI_ROOT_PATH, '{subfolder}', '{filename}.fastq.gz')
    conda:
        '../../environment/conda/conda_pbtools.yml'
    wildcard_constraints:
        filename = CONSTRAINT_BAM
    params:
        out_prefix = lambda wildcards, output: output.fastq.rsplit('.', 2)[0]
    resources:
        runtime_hrs = lambda wildcards, attempt: 11 * attempt
    shell:
        'bam2fastq -c 6 -o {params.out_prefix} {input.bam}'


rule compress_fastq:
    input:
        fastq = os.path.join(HIFI_ROOT_PATH, '{subfolder}', '{filename}.fastq')
    output:
        fastq = os.path.join(HIFI_ROOT_PATH, '{subfolder}', '{filename}.fastq.gz')
    conda:
        '../../environment/conda/conda_shelltools.yml'
    wildcard_constraints:
        filename = CONSTRAINT_FASTQ
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: 11 * attempt
    shell:
        'pigz --best -p {threads} -c {input.fastq} > {output.fastq}'