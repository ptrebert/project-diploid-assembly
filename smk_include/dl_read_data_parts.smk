
import collections

READ_SOURCE_URL = dict()
COMPRESSED_PARTS = []
UNCOMPRESSED_PARTS = []
PARTITIONED_SAMPLES = []
PARTITIONS_PER_SAMPLE = collections.defaultdict(list)

# needed below for merging
PARTS_DATA_PATH = os.path.join(os.getcwd(), 'input', 'read_data', 'parts')

for project, source_files in config['read_files'].items():
    for file_num, file_infos in source_files.items():
        if (bool(file_infos['skip'])):
            continue
        if file_num == 'fileN':
            # a chunked dataset is not
            # handled here
            continue

        local_name = file_infos['name']
        sample = local_name.rsplit('.', 2)[0]
        remote_url = file_infos['url']
        if remote_url.endswith('.gz'):
            COMPRESSED_PARTS.append(local_name)
        else:
            raise RuntimeError('Pipeline cannot handle remote uncompressed parts')
            #UNCOMPRESSED_PARTS.append(local_name)
        part_num = local_name.split('.')[-1]

        PARTITIONS_PER_SAMPLE[sample].append(part_num)

        assert local_name not in READ_SOURCE_URL
        READ_SOURCE_URL[local_name] = remote_url
        PARTITIONED_SAMPLES.append(sample)

PARTITIONED_SAMPLES = sorted(set(PARTITIONED_SAMPLES))

#rule prepare_partitioned_read_data:
#    input:
#        expand('input/read_data/diploid_assembly_input/{sample}.fastq.gz',
#                sample=PARTITIONED_SAMPLES)


rule download_compressed_read_parts:
    output:
        'input/read_data/parts/{part}.fastq.gz'
    log: 'log/download/read_parts/{part}.log'
    threads: 8
    run:
        exec = 'aria2c -s {threads} -x {threads}'
        exec += ' -o {output}'
        exec += ' {}'.format(READ_SOURCE_URL[wildcards.part])
        exec += ' &> {log}'
        shell(exec)


rule merge_partitioned_read_files:
    input:
        lambda wildcards: expand('input/read_data/parts/{sample}.part.{partnum}.fastq.gz',
                                 partnum=PARTITIONS_PER_SAMPLE[wildcards.sample], sample=wildcards.sample)
    output:
        'input/read_data/diploid_assembly_input/{sample}.fastq.gz'
    log: 'log/merge_parts/merge_{sample}.log'
    shell:
        "cat {input} > {output} 2> {log}"
