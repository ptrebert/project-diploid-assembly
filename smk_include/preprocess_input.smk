
include: 'handle_data_download.smk'

localrules: master_preprocess_input

rule master_preprocess_input:
    input:
        rules.master_handle_data_download.input,

        expand('input/fastq/complete/{sample}_1000.fastq.gz',
                sample=config['samples']),

        expand('input/fastq/complete/{sample}_sseq.fastq.gz',
                sample=config['strandseq_samples'])


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

    glob_pattern = os.path.join(checkpoint_dir, '{}_{{run_id}}.request'.format(sample))

    checkpoint_wildcards = glob_wildcards(glob_pattern)

    checkpoint_root = os.path.split(os.path.split(checkpoint_dir)[0])[0]
    fastq_input = os.path.join(checkpoint_root, '{individual}_{bioproject}', '{sample}_{run_id}.fastq.gz')

    fastq_files = expand(
        fastq_input,
        individual=individual,
        bioproject=bioproject,
        sample=sample,
        run_id=checkpoint_wildcards.run_id
        )

    return fastq_files


rule merge_strandseq_libraries:
    input:
        collect_strandseq_libraries
    output:
        'input/fastq/complete/{sample}_sseq.fastq.gz'
    log:
        'log/input/fastq/complete/{sample}_sseq.merge.log'
    run:
        with open(log[0], 'w') as logfile:
            _ = logfile.write('Start merging of {} files\n'.format(len(input)))
            try:
                with open(output[0], 'wb') as merge_file:
                    for fastq in input:
                        _ = logfile.write('Processing {}\n'.format(fastq))
                        _ = merge_file.write(open(fastq, 'rb').read())
                        _ = logfile.write('Completed: {}\n'.format(fastq))
            except Exception as err:
                _ = logfile.write('Error: {}'.format(str(err)))
                raise err
    # end of rule
