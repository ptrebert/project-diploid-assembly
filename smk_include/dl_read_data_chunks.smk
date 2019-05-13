
READ_SOURCE_URL = dict()
COMPRESSED_CHUNKS = []
UNCOMPRESSED_CHUNKS = []
CHUNKED_SAMPLES = []

# needed below for merging
CHUNK_DATA_PATH = os.path.join(os.getcwd(), 'input', 'read_data', 'chunks')

for project, source_files in config['read_files'].items():
    for file_num, file_infos in source_files.items():
        if (bool(file_infos['skip'])):
            continue
        if file_num != 'fileN':
            # not a dataset that comes in chunks
            # don't handle that here
            continue

        num_range = int(file_infos['range'])
        for n in range(1, num_range + 1):
            local_name = file_infos['name'].format(**{'FILENUM': str(n)})
            sample = local_name.rsplit('.', 2)[0]
            remote_url = file_infos['url'].format(**{'FILENUM': str(n)})
            if remote_url.endswith('.gz'):
                raise RuntimeError('Pipeline cannot handle remote compressed chunks')
                #COMPRESSED_CHUNKS.append(local_name)
            else:
                UNCOMPRESSED_CHUNKS.append(local_name)
            assert local_name not in READ_SOURCE_URL
            READ_SOURCE_URL[local_name] = remote_url
            CHUNKED_SAMPLES.append(sample)

CHUNKED_SAMPLES = sorted(set(CHUNKED_SAMPLES))

#rule prepare_chunked_read_data:
#    input:
#        expand('input/read_data/diploid_assembly_input/{sample}.fastq.gz',
#                sample=CHUNKED_SAMPLES)


rule download_uncompressed_reads_chunk:
    """
    Think twice before executing this rule 'somewhere'!
    May trigger the download of hundreds of small chunks
    in parallel and cause a s***load of I/O
    """
    output:
        'input/read_data/chunks/{chunk}.fastq.gz'
    log: 'log/download/read_chunks/{chunk}.log'
    threads: 1
    run:
        exec = 'wget --quiet -O - '
        exec += ' {}'.format(READ_SOURCE_URL[wildcards.chunk])
        exec += ' | gzip > {output}'
        exec += ' 2> {log}'
        shell(exec)


def collect_read_chunk_files(wildcards):

    subset = int(wildcards.subset)
    split = int(wildcards.split)

    filter_fun = lambda x: wildcards.sample in x and x.endswith('fastq.gz')
    input_files = []

    read_files = sorted(filter(filter_fun, os.listdir(CHUNK_DATA_PATH)), key=lambda x: int(x.split('.')[-3]))
    subset_size = int(len(read_files) // split)
    # this is left-exclusive as chunks are enumerated starting from 1
    lower_bound = (subset - 1) * subset_size
    upper_bound = subset * subset_size
    if subset == split:
        # correct for int taking only integral part (floor rounding)
        upper_bound += split
    for rf in read_files:
        # after validation, replace for loop with list slicing
        chunk_id = int(rf.split('.')[-3])
        if lower_bound < chunk_id <= upper_bound:
            input_files.append(os.path.join(CHUNK_DATA_PATH, rf))
    assert input_files, 'No read data part files collected for: {}'.format(wildcards)
    return input_files


rule merge_chunked_read_files:
    input:
        read_data = collect_read_chunk_files
    output:
        temp('input/read_data/merged_chunks/{sample}/{sample}.part.{subset}{split}.fastq.gz')
    log: 'log/merge_chunks/merge_{sample}.part.{subset}{split}.log'
    shell:
        "cat {input} > {output} 2> {log}"


rule merge_chunk_part_read_files:
    input:
        expand('input/read_data/merged_chunks/{{sample}}/{{sample}}.part.{subset}{split}.fastq.gz',
                subset=[1, 2, 3, 4, 5], split=[5])
    output:
        'input/read_data/diploid_assembly_input/{sample}.fastq.gz'
    log: 'log/merge_parts/merge_{sample}.log'
    shell:
        "cat {input} > {output} 2> {log}"
