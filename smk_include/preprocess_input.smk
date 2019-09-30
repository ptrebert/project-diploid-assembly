
include: 'handle_data_download.smk'

localrules: master_preprocess_input, ex_nihilo_input_injections

rule ex_nihilo_input_injections:
    output:
        protected(expand('input/fastq/complete/{individual}_hgsvc_pbsq2-ccs_1000.fastq.gz',
                         individual=['HG00733', 'HG00732', 'HG00731'])),
        protected('input/fastq/complete/HG002_pbio_pbsq2-ccs_1000.fastq.gz'),
        protected('input/fastq/complete/HG002_v19_pbsq2-ccs_1925.fastq.gz'),
        protected('input/fastq/complete/HG002_v20_pbsq2-ccs_1520.fastq.gz'),
        protected('input/bam/complete/HG002_pbio_pbsq2-clr_1000.pbn.bam'),


rule master_preprocess_input:
    input:
        rules.master_handle_data_download.input,
        rules.ex_nihilo_input_injections.output


def collect_fastq_input_parts(wildcards):

    base_path = 'input/fastq/partial/parts'

    sample = wildcards.sample

    all_part_files = os.listdir(base_path)
    relevant_parts = list(filter(lambda x: sample in x and x.endswith('fastq.gz'), all_part_files))

    relevant_parts = [os.path.join(base_path, f) for f in relevant_parts]

    if len(relevant_parts) == 0:
        raise RuntimeError('Empty input')

    return sorted(relevant_parts)


rule merge_fastq_input_parts:
    input:
        collect_fastq_input_parts
    output:
        'input/fastq/complete/{sample}_1000.fastq.gz'
    wildcard_constraints:
        sample = '(' + '|'.join(config['partial_samples']) + ')'
    log:
        'log/input/fastq/complete/{sample}_1000.merge.log'
    shell:
        'cat {input} > {output} 2> {log}'


def collect_strandseq_libraries(wildcards):

    sample = wildcards.sample
    bioproject = config['strandseq_to_bioproject'][sample]
    individual = sample.split('_')[0]

    checkpoint_dir = checkpoints.create_bioproject_download_requests.get(individual=individual, bioproject=bioproject).output[0]

    glob_pattern = os.path.join(checkpoint_dir, '{lib_id}.request')

    checkpoint_wildcards = glob_wildcards(glob_pattern)

    checkpoint_root = os.path.split(os.path.split(checkpoint_dir)[0])[0]
    fastq_input = os.path.join(checkpoint_root, '{individual}_{bioproject}', '{lib_id}.fastq.gz')

    fastq_files = expand(
        fastq_input,
        individual=individual,
        bioproject=bioproject,
        lib_id=checkpoint_wildcards.lib_id
        )

    return fastq_files


rule merge_strandseq_libraries:
    """
    TODO
    This creates a mock file for
    the time being to keep the pipeline
    in a coherent state while the strand-seq
    steps are being implemented
    """
    input:
        collect_strandseq_libraries
    output:
        'input/fastq/complete/{sample}_sseq.fastq.gz'
    message: 'Creating strand-seq mock file'
    run:
        with open(output[0], 'w') as dump:
            for fastq in sorted(input):
                dump.write(fastq + '\n')


def collect_pacbio_bam_input_parts(wildcards):

    base_path = 'input/bam/partial/parts'

    sample = wildcards.sample

    all_part_files = os.listdir(base_path)
    relevant_parts = list(filter(lambda x: sample in x and x.endswith('pbn.bam'), all_part_files))

    relevant_parts = [os.path.join(base_path, f) for f in relevant_parts]

    if len(relevant_parts) == 0:
        raise RuntimeError('Empty input')

    return sorted(relevant_parts)


rule merge_pacbio_native_bams:
    input:
        parts = collect_pacbio_bam_input_parts
    output:
        'input/bam/complete/{sample}_1000.pbn.bam'
    log:
        'log/input/bam/complete/{sample}_1000.mrg.log'
    benchmark:
        'run/input/bam/complete/{sample}_1000.mrg.rsrc'
    params:
        input_list = lambda wildcards, input: ' -in ' + ' -in '.join(input.parts)
    shell:
        'bamtools merge {params.input_list} -out {output}'