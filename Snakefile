
import collections as col

wildcard_constraints:
    sample="HG[0-9]{3,5}"

REFKEY = config['use_reference_genome']
REFNAME = config['ref_genome'][REFKEY]['name']
REFCHROM = config['chromosomes'][REFKEY]

# Collect ALL
READ_FILES = []
READ_SOURCE_URL = dict()

DOWNLOAD_SORT = col.defaultdict(list)

SAMPLES = []
PROJECTS = []

def get_file_sort_key(filename, file_url):

    status = 1 if file_url.endswith('.gz') else 0
    if 'part' in filename:
        key = (0, 1, 0, status)
    elif 'chunk' in filename:
        key = (1, 0, 0, status)
    elif 'cmp' in filename:
        key = (0, 0, 1, status)
    else:
        raise ValueError('Could not determine sort key for file: {}'.format(filename))
    return key


for project, source_files in config['read_files'].items():
    for file_num, file_infos in source_files.items():
        if (bool(file_infos['skip'])):
            continue
        if file_num != 'fileN':
            local_name = file_infos['name']
            sample, project, _ = local_name.split('.', 2)

            sort_key = get_file_sort_key(local_name, file_infos['url'])
            DOWNLOAD_SORT[sort_key].append(local_name)

            SAMPLES.append(sample)
            PROJECTS.append(project)
            READ_FILES.append(local_name)
            READ_SOURCE_URL[local_name] = file_infos['url']
        else:
            num_range = int(file_infos['range'])
            for n in range(1, num_range + 1):
                local_name = file_infos['name'].format(**{'FILENUM': str(n)})
                sample, project, _ = local_name.split('.', 2)
                sort_key = get_file_sort_key(local_name, file_infos['url'])
                remote_url = file_infos['url'].format(**{'FILENUM': str(n)})

                sort_key = get_file_sort_key(local_name, file_infos['url'])
                DOWNLOAD_SORT[sort_key].append(local_name)

                SAMPLES.append(sample)
                PROJECTS.append(project)
                READ_FILES.append(local_name)
                READ_SOURCE_URL[local_name] = remote_url

SAMPLES = sorted(set(SAMPLES))
PROJECTS = sorted(set(PROJECTS))


rule master:
    input:
        'references/' + REFNAME,
        # Output files: download tasks
        expand('input/read_data/part/{gzfile}.fastq.gz',
                gzfile=DOWNLOAD_SORT[(0, 1, 0, 1)]),

        expand('input/read_data/chunk/{rawfile}.fastq.gz',
                rawfile=DOWNLOAD_SORT[(1, 0, 0, 0)]),

        expand('input/read_data/cmp/{rawfile}.fastq.gz',
                rawfile=DOWNLOAD_SORT[(0, 0, 1, 0)]),

        # Mapping of chunks fails w/o error message - merge chunks
        # before mapping...
        expand('input/read_data/mrg/{sample}.{project}.ont-ul.part.{subset}{split}.fastq.gz',
                subset=[1, 2, 3, 4, 5], split=[5],
                sample=['HG002'], project=['ucsc1']),

        # Output files: read mapping against reference
        expand('output/alignments/{sample}.{project}.ont-ul.cram',
                sample=['HG002', 'HG00733'], project=['pangen']),
        expand('output/alignments/{sample}.{project}.ont-ul.cram',
                sample=['HG002'], project=['giab']),
        expand('output/alignments/{sample}.{project}.ont-ul.cram',
                sample=['HG002'], project=['ucsc1']),

        expand('output/sorted_aln/{sample}.{project}.ont-ul.sorted.cram',
                sample=['HG002', 'HG00733'], project=['pangen']),
        expand('output/sorted_aln/{sample}.{project}.ont-ul.sorted.cram',
                sample=['HG002'], project=['giab']),
        expand('output/sorted_aln/{sample}.{project}.ont-ul.sorted.cram',
                sample=['HG002'], project=['ucsc1']),

        # Phase variants
        expand('output/phased_variants/{sample}.{project}.{chromosome}.vcf.gz',
                sample=['HG00733'], project=['pangen'],
                chromosome=REFCHROM),

        # Merge phased variants
        expand('output/merged_phased_variants/{sample}.{project}.vcf.gz',
                sample=['HG00733'], project=['pangen']),

        # Perform haplo-tagging
        expand('output/tagged_aln/{sample}.{project}.ont-ul.sorted.tagged.cram',
                sample=['HG00733'], project=['pangen']),
        expand('output/tagged_aln/{sample}.{project}.ont-ul.haplotag.tsv.gz',
                sample=['HG00733'], project=['pangen'])

#        expand('phased-snvs/pb/HG00733.{chromosome}.vcf.gz', chromosome=chromosomes),
#        expand('phased-snvs/pb_ss/HG00733.{chromosome}.vcf.gz', chromosome=chromosomes),
#        'pacbio/cram/HG00733.sorted.tagged.cram',
#        'pacbio/fastq/HG00733.h1-only.fastq.gz',
#        'pacbio/fastq/HG00733.h2-only.fastq.gz',
#        'pacbio/fastq/HG00733.h1-plus-untagged.fastq.gz',
#        'pacbio/fastq/HG00733.h2-plus-untagged.fastq.gz',
    message: 'Executing ALL'


onsuccess:
    body_text = "Nothing to report\n{log}"
    if config['notify']:
        shell('mail -s "[Snakemake] diploid genome assembly - success" {} <<< "{}"'.format(config['notify_email'], body_text))

onerror:
    if config['notify']:
        shell('mail -s "[Snakemake] diploid genome assembly - ERRROR" {} < {{log}}'.format(config['notify_email']))


rule download_reference_genome:
    output:
        'references/' + REFNAME,
    log: 'log/DL_{}.log'.format(REFNAME)
    threads: 4
    params:
        seq_url = config['ref_genome'][REFKEY]['url']
    run:
        exec = 'aria2c -s {threads} -x {threads}'
        exec += ' -o {output}'
        exec += ' {params.seq_url}'
        exec += ' &> {log}'
        shell(exec)


rule download_big_gzipped_reads:
    output:
        expand('input/read_data/{size}/{{gzfile}}.fastq.gz',
                size=['part'])
    log: 'log/download/{gzfile}.log'
    threads: 8
    run:
        exec = 'aria2c -s {threads} -x {threads}'
        exec += ' -o {output}'
        exec += ' {}'.format(READ_SOURCE_URL[wildcards.gzfile])
        exec += ' &> {log}'
        shell(exec)


rule download_uncompressed_reads_chunk:
    """
    Think twice before executing this rule 'somewhere'!
    May trigger the download of hundreds of small chunks
    in parallel and cause a s***load of I/O
    """
    output:
        expand('input/read_data/{size}/{{rawfile}}.fastq.gz',
                size=['chunk']),
    log: 'log/download/{rawfile}.log'
    threads: 1
    run:
        exec = 'wget --quiet -O - '
        exec += ' {}'.format(READ_SOURCE_URL[wildcards.rawfile])
        exec += ' | gzip > {output}'
        exec += ' 2> {log}'
        shell(exec)


rule download_uncompressed_reads_complete:
    output:
        expand('input/read_data/{size}/{{rawfile}}.fastq.gz',
                size=['cmp'])
    log: 'log/download/{rawfile}.log'
    threads: 1
    run:
        exec = 'wget --quiet -O - '
        exec += ' {}'.format(READ_SOURCE_URL[wildcards.rawfile])
        exec += ' | gzip > {output}'
        exec += ' 2> {log}'
        shell(exec)

def collect_read_files_merging(wildcards):

    subset = int(wildcards.subset)
    split = int(wildcards.split)

    filter_fun = lambda x: wildcards.sample in x and wildcards.project in x
    data_path = os.path.join(os.getcwd(), 'input', 'read_data')
    input_files = []
    for root, dirs, files in os.walk(data_path):
        if any([root.endswith(s) for s in ['part', 'chunk']]):
            read_files = sorted(filter(filter_fun, files), key=lambda x: int(x.split('.')[-3]))
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
                    input_files.append(os.path.join(root, rf))
    assert input_files, 'No read data part files collected for: {}'.format(wildcards)
    return input_files


rule merge_chunked_read_files:
    input:
        read_data = collect_read_files_merging
    output:
        'input/read_data/mrg/{sample}.{project}.ont-ul.part.{subset}{split}.fastq.gz',
    log: 'log/merge_{sample}.{project}.{subset}{split}.log'
    shell:
        "cat {input} > {output} 2> {log}"


def collect_read_files_alignment(wildcards):

    filter_fun = lambda x: wildcards.sample in x and wildcards.project in x
    data_path = os.path.join(os.getcwd(), 'input', 'read_data')
    input_files = []
    for root, dirs, files in os.walk(data_path):
        if any([root.endswith(s) for s in ['part', 'cmp', 'mrg']]):
            read_files = filter(filter_fun, files)
            read_files = [os.path.join(root, f) for f in sorted(read_files)]
            input_files.extend(read_files)
    #assert input_files, 'No read data part files collected for: {} / {}'.format(wildcards.sample, wildcards.project)
    return input_files


rule map_reads:
    input:
        reference = 'references/' + REFNAME,
        read_data = collect_read_files_alignment
    output:
        'output/alignments/{sample}.{project}.ont-ul.cram'
    log:
        minimap = 'log/{sample}.{project}.minimap.log',
        samtools = 'log/{sample}.{project}.samtools.log'
    params:
        cache_path = os.path.join(os.getcwd(), 'references', 'cache')
    threads: 48
    run:
        exec = 'minimap2 -t {threads}'
        exec += ' -R \'@RG\\tID:1\\tSM:{wildcards.sample}\''
        exec += ' -a -x map-ont'
        exec += ' {input.reference}'
        exec += ' {input.read_data} 2> {log.minimap}'
        exec += ' | samtools view'
        exec += ' -T {input.reference}'
        exec += ' -C -'
        exec += ' > {output}'
        exec += ' 2> {log.samtools}'
        shell(exec)


# sorting takes ~1 hour / 50 GiB with 20 cores and 5GiB / core
# ATTENTION: current settings result in 20 * 20 GiB RAM => 400 GiB required!
rule sort_cram_alignments:
    input:
        'output/alignments/{sample}.{project}.ont-ul.cram'
    output:
        'output/sorted_aln/{sample}.{project}.ont-ul.sorted.cram'
    log: 'log/aln_sort/{sample}.{project}.sort.log'
    threads: 20
    params:
        mem_limit = '20G',
        cache_path = os.path.join(os.getcwd(), 'references', 'cache')
    run:
        exec = 'REF_CACHE={params.cache_path} '
        exec += 'samtools sort -o {output} -m {params.mem_limit}'
        exec += ' -@ {threads} {input} &> {log}'
        shell(exec)


rule index_cram_alignments:
    input:
        cram='{cramfile}.sorted.cram'
    output:
        crai='{cramfile}.sorted.cram.crai'
    threads: 8
    shell:
        "samtools index -@ {threads} {input.cram}"


rule phase_variants:
    input:
        vcf = 'output/variant_calls/{sample}/{sample}.{chromosome}.vcf.gz',
        tbi = 'output/variant_calls/{sample}/{sample}.{chromosome}.vcf.gz.tbi',
        cram = 'output/sorted_aln/{sample}.{project}.ont-ul.sorted.cram',
        crai = 'output/sorted_aln/{sample}.{project}.ont-ul.sorted.cram.crai',
        ref = 'references/' + REFNAME,
    output:
        vcf = 'output/phased_variants/{sample}.{project}.{chromosome}.vcf.gz',
    log: 'log/phase_variants/{sample}.{project}.{chromosome}.vcf.log'
    run:
        exec = 'whatshap phase --reference {input.ref}'
        exec += ' --chromosome {wildcards.chromosome}'
        exec += ' --output {output.vcf}'
        exec += ' {input.vcf} {input.cram}'
        exec += ' &> {log}'
        shell(exec)


rule merge_phased_variants:
    input:
        expand('output/phased_variants/{{sample}}.{{project}}.{chromosome}.vcf.gz',
                chromosome=REFCHROM)
    output:
        vcf = 'output/merged_phased_variants/{sample}.{project}.vcf.gz',
        tbi = 'output/merged_phased_variants/{sample}.{project}.vcf.gz.tbi'
    log: 'log/merge_phased/{sample}.{project}.log'
    run:
        exec = 'bcftools concat -o {output.vcf} -O z {input} &> {log}'
        exec += ' && '
        exec += 'bcftools index --tbi --output-file {output.tbi} {output.vcf} &>> {log}'
        shell(exec)


rule haplotag_reads:
    input:
        vcf = 'output/merged_phased_variants/{sample}.{project}.vcf.gz',
        tbi = 'output/merged_phased_variants/{sample}.{project}.vcf.gz.tbi',
        cram = 'output/sorted_aln/{sample}.{project}.ont-ul.sorted.cram',
        crai = 'output/sorted_aln/{sample}.{project}.ont-ul.sorted.cram.crai',
        ref = 'references/' + REFNAME,
    output:
        cram = 'output/tagged_aln/{sample}.{project}.ont-ul.sorted.tagged.cram',
        taglist = 'output/tagged_aln/{sample}.{project}.ont-ul.haplotag.tsv.gz'
    log: 'log/tagged_aln/{sample}.{project}.haplotag.log'
    conda:
        "environment/conda/wh_split.yml"
    shell:
        "whatshap haplotag --output {output.cram} --reference {input.ref} --output-haplotag-list {output.taglist} {input.vcf} {input.cram} &> {log}"
    #run:
    #    exec = 'whatshap haplotag --output {output.cram}'
    #    exec += ' --reference {input.ref} --output-haplotag-list {output.taglist}'
    #    exec += ' {input.vcf} {input.cram}'
    #    exec += ' &> {log}'
    #    shell(exec)
