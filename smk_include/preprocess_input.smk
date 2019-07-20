
include: 'handle_data_download.smk'

localrules: master_preprocess_input

rule master_preprocess_input:
    input:
        rules.master_handle_data_download.input,

        expand('input/fastq/complete/{sample}_1000.fastq.gz',
                sample=config['samples'])


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
    log:
        'log/input/fastq/complete/{sample}_1000.merge.log'
    shell:
        'cat {input} > {output} 2> {log}'