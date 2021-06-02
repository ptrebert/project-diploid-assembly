localrules: master

# config values data
SAMPLES = config['samples']
RUN_SYSTEM = config['system']

SHORT_READS_PATH = '' # config['short_reads_path'][SAMPLE][RUN_SYSTEM]

REFERENCE_FOLDER = config['reference_folder'][RUN_SYSTEM]
REFERENCE_ASSEMBLIES = config['reference_assemblies']

# config values system
RUNTIME_ES = config['runtime']['ExtractSequence'][RUN_SYSTEM]

MEMORY_HIFIASM = config['memory']['HifiAsm']

# config values parameters
MBG_KMER_SIZE = config['params']['MBG']['kmer']
MBG_WINDOW_SIZE = config['params']['MBG']['window']

GA_MIN_SCORE = config['params']['GraphAligner']['minscore']

# size fraction zero = dump all reads
READ_SIZE_FRACTIONS = config['params']['SizeFractions']

MAIN_CHROMOSOMES = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']

ALL_ONT_FILES = {
    'HG002': [
        'HG002_giab_ULfastqs_guppy324',
        'HG002_ONT_PAD64459_Guppy32',
    ],

    'HG00733': [
        'HG00733_EEE_ONT',
        'HG00733_1_Guppy_4.2.2_prom',
        'HG00733_2_Guppy_4.2.2_prom',
        'HG00733_3_Guppy_4.2.2_prom'
    ],
}

def find_script_path(script_name, subfolder=''):
    """
    Find full path to script to be executed. Function exists
    to avoid config parameter "script_dir"

    :param script_name:
    :param subfolder:
    :return:
    """
    import os

    current_root = workflow.basedir
    last_root = ''

    script_path = None

    for _ in range(workflow.basedir.count('/')):
        if last_root.endswith('project-diploid-assembly'):
            raise RuntimeError('Leaving project directory tree (next: {}). '
                               'Cannot find script {} (subfolder: {}).'.format(current_root, script_name, subfolder))
        check_path = os.path.join(current_root, 'scripts', subfolder).rstrip('/')  # if subfolder is empty string
        if os.path.isdir(check_path):
            check_script = os.path.join(check_path, script_name)
            if os.path.isfile(check_script):
                script_path = check_script
                break
        last_root = current_root
        current_root = os.path.split(current_root)[0]

    if script_path is None:
        raise RuntimeError('Could not find script {} (subfolder {}). '
                           'Started at path: {}'.format(script_name, subfolder, workflow.basedir))
    return script_path


rule run_all:
    input:
        expand(
            'output/read_info/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.h5',
            filename=ALL_ONT_FILES['HG00733'],
            kmer=MBG_KMER_SIZE * 4,
            window=MBG_WINDOW_SIZE * 4,
            minscore=[GA_MIN_SCORE] * 4 * len(MBG_KMER_SIZE)
        ),
        expand(
            'output/read_info/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.h5',
            filename=ALL_ONT_FILES['HG002'],
            kmer=MBG_KMER_SIZE * 2,
            window=MBG_WINDOW_SIZE * 2,
            minscore=[GA_MIN_SCORE] * 2 * len(MBG_KMER_SIZE)
        )


def select_raw_ont_files(wildcards):

    if wildcards.filename == 'HG00733_EEE_ONT':
        input_folder = config['ont_raw_folder']['HG00733'][RUN_SYSTEM]
        input_files = config['ont_split_files']
        ont_files = [os.path.join(input_folder, f + '.fastq.gz') for f in input_files]
    elif 'HG00733' in wildcards.filename:
        input_folder = config['ont_raw_folder']['HG00733'][RUN_SYSTEM]
        ont_files = [os.path.join(input_folder, wildcards.filename + '.fastq.gz')]
    elif 'HG002' in wildcards.filename:
        input_folder = config['ont_raw_folder']['HG002'][RUN_SYSTEM]
        ont_files = [os.path.join(input_folder, wildcards.filename + '.fastq.gz')]
    else:
        raise ValueError('Cannot handle ONT input request: {}'.format(str(wildcards)))
    return ont_files


rule clean_hpg_ont:
    input:
        fastq = select_raw_ont_files
    output:
        fasta = 'input/ont/{filename}.fa.gz'
    benchmark:
        'rsrc/input/ont/{filename}.clean.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt ** attempt
    params:
        pigz_cpu = lambda wildcards: int(config['num_cpu_medium']) // 2
    shell:
        'pigz -p {params.pigz_cpu} -c -d {input.fastq} | '
        'seqtk seq -S -A -C {input.fastq} | '
        'pigz -p {params.pigz_cpu} > {output.fasta}'


def set_mbg_memory(wildcards, attempt):

    if int(wildcards.kmer) < 1000:
        return config['memory']['MBG_high'][wildcards.sample] * attempt
    else:
        return config['memory']['MBG_low'][wildcards.sample] * attempt


def select_hifi_input(wildcards):

    return config['hifi_reads_path'][wildcards.sample][RUN_SYSTEM]


rule build_hifi_read_graph:
    """
    old call for MBG v1.03:
    MBG -i {input} -o {output} -t {threads} --blunt -k {wildcards.kmer} -w {wildcards.window} &> {log}

    MBG v1.04+ contains fix for overlap-longer-than-nodes error
    """
    input:
        select_hifi_input
    output:
        'output/mbg_hifi/{sample}_HiFi.mbg-k{kmer}-w{window}.gfa'
    log:
        'log/output/mbg_hifi/{sample}_HiFi.mbg-k{kmer}-w{window}.log'
    benchmark:
        'rsrc/output/mbg_hifi/{sample}_HiFi.mbg-k{kmer}-w{window}.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = set_mbg_memory,
        runtime_hrs = lambda wildcards, attempt: attempt + attempt * attempt
    shell:
        'MBG -i {input} -o {output} -t {threads} -k {wildcards.kmer} -w {wildcards.window} &> {log}'


rule clean_mbg_graph:
    """
    https://github.com/maickrau/MBG/issues/1

    ### I've used vg: vg view -Fv graph.gfa | vg mod -n -U 100 - | vg view - > blunt-graph.gfa
    """
    input:
        'output/mbg_hifi/{sample}_HiFi.mbg-k{kmer}-w{window}.gfa'
    output:
        'output/mbg_hifi_clean/{sample}_HiFi.mbg-k{kmer}-w{window}.clean.gfa'
    log:
        'log/output/mbg_hifi_clean/{sample}_HiFi.mbg-k{kmer}-w{window}.clean.log'
    benchmark:
        'rsrc/output/mbg_hifi_clean/{sample}_HiFi.mbg-k{kmer}-w{window}.clean.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt ** attempt
    shell:
        'vg view -Fv {input} | vg mod -n -U 100 - 2> {log} | vg view - > {output}'


rule extract_roi_sequences:
    input:
        ref_fasta = os.path.join(REFERENCE_FOLDER, '{reference}.fasta'),
        bed = os.path.join(REFERENCE_FOLDER, '{reference}.{subset}.bed'),
    output:
        fasta = os.path.join(REFERENCE_FOLDER, '{reference}.{subset}.fasta'),
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'bedtools getfasta -nameOnly -fo {output.fasta} -fi {input.ref_fasta} -bed {input.bed}'


rule align_roi_sequences_to_mbg:
    input:
        graph = 'output/mbg_hifi/{sample}_HiFi.mbg-k{kmer}-w{window}.gfa',
        reads = os.path.join(REFERENCE_FOLDER, '{reference}.{subset}.fasta'),
    output:
        gaf = 'output/roi_aln/{reference}.{subset}_MAP-TO_{sample}_HiFi.mbg-k{kmer}-w{window}.gaf',
    log:
        'log/output/roi_aln/{reference}.{subset}_MAP-TO_{sample}_HiFi.mbg-k{kmer}-w{window}.ga.log',
    benchmark:
        'rsrc/output/roi_aln/{reference}.{subset}_MAP-TO_{sample}_HiFi.mbg-k{kmer}-w{window}.ga.rsrc',
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 65536 * attempt,
        runtime_hrs = lambda wildcards, attempt: 12 * attempt
    shell:
        'GraphAligner -t {threads} -g {input.graph} -f {input.reads} '
            '-x dbg --min-alignment-score 500 --multimap-score-fraction 1 '
            '-a {output.gaf} &> {log}'


rule strip_sequences_from_graph:
    input:
        gfa = 'output/mbg_hifi/{sample}_HiFi.mbg-k{kmer}-w{window}.gfa'
    output:
        gfa = 'output/mbg_hifi/{sample}_HiFi.mbg-k{kmer}-w{window}.noseq.gfa'
    benchmark:
        'rsrc/output/mbg_hifi/{sample}_HiFi.mbg-k{kmer}-w{window}.noseq.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'gfatools view -S {input.gfa} > {output.gfa}'


def select_mbg_graph(wildcards):

    path = 'output/mbg_hifi/{sample}_HiFi.mbg-k{kmer}-w{window}.gfa'
    formatter = dict(wildcards)

    if 'HG002' in wildcards.filename or 'NA24385' in wildcards.filename:
        formatter['sample'] = 'HG002'
    elif 'HG00733' in wildcards.filename:
        formatter['sample'] = 'HG00733'
    else:
        raise ValueError('Cannot handle wildcards: {}'.format(formatter))

    return path.format(**formatter)


rule ont_error_correction:
    """
    """
    input:
        graph = select_mbg_graph,
        reads = 'input/ont/{filename}.fa.gz',
    output:
        gaf = 'output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.gaf',
        ec_reads_clip = 'output/ont_ec/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.clip-ec.fa.gz',
        ec_reads_raw = 'output/ont_ec/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.raw-ec.fa.gz',  # potentially useful for scaffolding
    log:
        'log/output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.ga.log'
    benchmark:
        'rsrc/output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.ga.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 65536 * attempt,
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    shell:
        'GraphAligner -g {input.graph} -f {input.reads} '
            '-x dbg -t {threads} '
            '--min-alignment-score {wildcards.minscore} --multimap-score-fraction 1 '
            '--corrected-clipped-out {output.ec_reads_clip} --corrected-out {output.ec_reads_raw}'
            '-a {output.gaf} &> {log}'


rule get_sequence_stats:
    input:
        '{filepath}/{filename}.fa.gz'
    output:
        'output/seq_stats/{filepath}/{filename}.stats.tsv.gz'
    benchmark:
        'rsrc/output/seq_stats/{filepath}/{filename}.seqtk.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        pigz_cpu = lambda wildcards: int(config['num_cpu_low']) // 2
    shell:
        'pigz -p {params.pigz_cpu} -d -c {input} | '
        'seqtk comp | '
        'pigz -p {params.pigz_cpu} --best > {output}'


rule compute_stats_input_reads:
    input:
        tsv = lambda wildcards: expand(
            'output/seq_stats/input/ont/{filename}.stats.tsv.gz',
            filename=ALL_ONT_FILES[wildcards.sample]
        )
    output:
        summary = 'output/seq_stats/summary/{sample}_ONT_input.stats.tsv',
        dump = 'output/seq_stats/summary/{sample}_ONT_input.stats.pck'
    log:
        'log/output/seq_stats/summary/{sample}_ONT_input.stats.log'
    benchmark:
        'rsrc/output/seq_stats/summary/{sample}_ONT_input.stats.rsrc'
    conda:
        '../../environment/conda/conda_pyscript.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
        '{params.script_exec} --debug --input-files {input.tsv} '
        '--output {output.dump} --summary-output {output.summary} '
        '--num-cpu {threads} --genome-size 3100000000 &> {log}'


rule extract_read_info_from_gaf:
    input:
        gaf = 'output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.gaf'
    output:
        hdf = 'output/read_info/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.h5',
        listings = expand(
            'output/read_info/{{filename}}_MAP-TO_mbg-k{{kmer}}-w{{window}}.ms{{minscore}}.geq{size_fraction}.{ext}',
            size_fraction=[0] + READ_SIZE_FRACTIONS,
            ext=['path-nodes.txt', 'readec-path.tsv', 'read-ec.txt']
            )
    log:
        'log/output/read_info/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.gaf-filter.log',
    benchmark:
        'rsrc/output/read_info/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.gaf-filter.rsrc',
    conda:
        '../../environment/conda/conda_pyscript.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: max(0, attempt-1) * attempt
    params:
        script_exec = lambda wildcards: find_script_path('filter_gaf.py'),
        size_fractions = lambda wildcards: ' '.join(list(map(str, READ_SIZE_FRACTIONS)))
    shell:
        '{params.script_exec} --debug --input {input.gaf} '
        '--output {output.hdf} --size-fractions {params.size_fractions} &> {log}'


rule deduplicate_ec_reads:
    input:
        listings = expand(
            'output/read_info/{{filename}}_MAP-TO_mbg-k{kmer}-w{window}.ms{{minscore}}.geq{{size_fraction}}.read-ec.txt',
            zip,
            kmer=MBG_KMER_SIZE,
            window=MBG_WINDOW_SIZE,
            )
    output:
        listings = expand(
            'output/read_info/{{filename}}_MAP-TO_mbg-k{kmer}-w{window}.ms{{minscore}}.geq{{size_fraction}}.read-ec.dedup.txt',
            zip,
            kmer=MBG_KMER_SIZE,
            window=MBG_WINDOW_SIZE,
            )
    benchmark:
        'rsrc/output/read_info/{filename}_MAP-TO_mbg-kALL-wALL.ms{minscore}.geq{size_fraction}.dedup.rsrc',
    run:
        import collections as col
        import random as rand
        import io as io

        assert isinstance(input.listings, list), 'Expected more than one input file: {}'.format(input)

        read_to_readset = col.defaultdict(list)
        output_buffers = dict()

        for read_ec_file in input.listings:
            output_file = read_ec_file.rsplit('.', 1)[0] + '.dedup.txt'
            output_buffers[output_file] = io.StringIO()
            with open(read_ec_file, 'r') as listing:
                for line in listing:
                    read_to_readset[line.strip()].append(output_file)
        
        for read, readset in read_to_readset.items():
            select_buffer = output_buffers[rand.choice(readset)]
            _ = select_buffer.write(read + '\n')

        for output_file, output_buffer in output_buffers.items():
            with open(output_file, 'w') as dump:
                _ = dump.write(output_buffer.getvalue())
        

rule extract_ec_reads_by_size:
    input:
        ec_reads = 'output/ont_ec/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.clip-ec.fa.gz',
        read_list = 'output/read_info/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.geq{size_fraction}.read-ec.dedup.txt'
    output:
        'output/ont_ec_subsets/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.clip-ec.geq{size_fraction}.fa.gz'
    benchmark:
        'rsrc/output/ont_ec_subsets/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.clip-ec.geq{size_fraction}.seqtk.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048,
        runtime_hrs = lambda wildcards, attempt: RUNTIME_ES * attempt
    params:
        pigz_cpu = lambda wildcards: int(config['num_cpu_low']) // 2
    shell:
        'pigz -p {params.pigz_cpu} -d -c {input.ec_reads} | '
        'seqtk subseq /dev/stdin {input.read_list} | '
        'pigz -p {params.pigz_cpu} --best > {output}'


rule compute_stats_corrected_reads:
    input:
        tsv = lambda wildcards: expand(
            'output/seq_stats/output/ont_ec_subsets/{filename}_MAP-TO_mbg-k{{kmer}}-w{{window}}.ms{{minscore}}.clip-ec.geq0.stats.tsv.gz',
            filename=ALL_ONT_FILES[wildcards.sample]
            )
    output:
        summary = 'output/seq_stats/summary/{sample}_ONT_clip-ec.k{kmer}-w{window}.ms{minscore}.stats.tsv',
        dump = 'output/seq_stats/summary/{sample}_ONT_clip-ec.k{kmer}-w{window}.ms{minscore}.stats.pck'
    log:
        'log/output/seq_stats/summary/{sample}_ONT_clip-ec.k{kmer}-w{window}.ms{minscore}.stats.log'
    benchmark:
        'rsrc/output/seq_stats/summary/{sample}_ONT_clip-ec.k{kmer}-w{window}.ms{minscore}.stats.rsrc'
    conda:
        '../../environment/conda/conda_pyscript.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
        '{params.script_exec} --debug --input-files {input.tsv} '
        '--output {output.dump} --summary-output {output.summary} '
        '--num-cpu {threads} --genome-size 3100000000 &> {log}'


def set_hifiasm_memory(wildcards, attempt, run_system):
    if attempt > 1:
        raise RuntimeError('hifiasm assembly failed first attempt: {}'.format(wildcards))
    if run_system == 'valet':
        return MEMORY_HIFIASM
    else:
        return MEMORY_HIFIASM * 2  # ~3 TB


# rule hifiasm_hifi_ontec_assembly:
#     """
#     """
#     input:
#         hifi_reads = select_hifi_input,
#         ontec_reads = expand(
#             'output/ont_ec_subsets/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.clip-ec.geq{{size_fraction}}.fa.gz',
#             zip,
#             filename=ONT_ALL_FILES * len(MBG_KMER_SIZE),
#             kmer=MBG_KMER_SIZE * len(ONT_ALL_FILES),
#             window=MBG_WINDOW_SIZE * len(ONT_ALL_FILES),
#             minscore=[GA_MIN_SCORE] * len(MBG_KMER_SIZE) * len(ONT_ALL_FILES)
#         )
#     output:
#         primary_unitigs = 'output/assembly/layout/{}_hifi-ontec_geq{{size_fraction}}/{}_hifi-ontec_geq{{size_fraction}}.p_utg.gfa'.format(SAMPLE, SAMPLE),
#         primary_contigs = 'output/assembly/layout/{}_hifi-ontec_geq{{size_fraction}}/{}_hifi-ontec_geq{{size_fraction}}.p_ctg.gfa'.format(SAMPLE, SAMPLE),
#         raw_unitigs = 'output/assembly/layout/{}_hifi-ontec_geq{{size_fraction}}/{}_hifi-ontec_geq{{size_fraction}}.r_utg.gfa'.format(SAMPLE, SAMPLE),
#         discard = multiext(
#             'output/assembly/layout/{}_hifi-ontec_geq{{size_fraction}}/{}_hifi-ontec_geq{{size_fraction}}'.format(SAMPLE, SAMPLE),
#             '.a_ctg.gfa', '.a_ctg.noseq.gfa',
#             '.ec.bin', '.ovlp.reverse.bin', '.ovlp.source.bin',
#             '.p_ctg.noseq.gfa', '.p_utg.noseq.gfa',
#             '.r_utg.noseq.gfa'
#         )
#     log:
#         hifiasm = 'log/output/assembly/layout/{}_hifi-ontec_{{size_fraction}}.hifiasm.log'.format(SAMPLE)
#     benchmark:
#         'rsrc/output/assembly/layout/{}_hifi-ontec_{{size_fraction}}.hifiasm'.format(SAMPLE) + '.t{}.rsrc'.format(config['num_cpu_max']) 
#     conda:
#         '../environment/conda/conda_biotools.yml'       
#     threads: config['num_cpu_max']
#     resources:
#         mem_total_mb = lambda wildcards, attempt: set_hifiasm_memory(wildcards, attempt, RUN_SYSTEM),
#         runtime_hrs = lambda wildcards, attempt: 52 + 52 * attempt
#     params:
#         prefix = lambda wildcards, output: output.primary_contigs.rsplit('.', 2)[0],
#     shell:
#         'hifiasm -o {params.prefix} -t {threads} --primary {input.hifi_reads} {input.ontec_reads} &> {log.hifiasm}'


rule count_reference_kmers:
    """
    k-mer sizes according to winnowmap github:
    15: read-to-ref mapping
    19: assm-to-ref mapping
    """
    input:
        fasta = os.path.join(REFERENCE_FOLDER, '{reference}.fasta'),
    output:
        kmer_db = directory(os.path.join(REFERENCE_FOLDER, '{reference}.k{kmer_size}.db/')),
        rep_kmer = os.path.join(REFERENCE_FOLDER, '{reference}.k{kmer_size}.rep-grt09998.txt')
    benchmark:
        'rsrc/output/kmer/{reference}.k{kmer_size}.count-dump.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(32768 * attempt / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        mem_total_gb = lambda wildcards, attempt: 32 * attempt
    shell:
        'meryl count k={wildcards.kmer_size} threads={threads} memory={resources.mem_total_gb} output {output.kmer_db} {input.fasta} && '
        'meryl print greater-than distinct=0.9998 {output.kmer_db} > {output.rep_kmer}'


# rule align_hifi_input_reads:
#     input:
#         reads = select_hifi_input,
#         fasta = os.path.join(REFERENCE_FOLDER, '{reference}.fasta'),
#         rep_kmer = os.path.join(REFERENCE_FOLDER, '{reference}.k15.rep-grt09998.txt')
#     output:
#         bam = 'output/read_align/{}_HiFi_input_MAP-TO_{{reference}}.k15.bam'.format(SAMPLE)
#     benchmark:
#         'rsrc/output/read_align/{}_HiFi_input_MAP-TO_{{reference}}.k15.wmap.rsrc'.format(SAMPLE)
#     log:
#         'log/output/read_align/{}_HiFi_input_MAP-TO_{{reference}}.k15.wmap.log'.format(SAMPLE)
#     conda:
#         '../../environment/conda/conda_biotools.yml'
#     threads: config['num_cpu_high']
#     resources:
#         mem_total_mb = lambda wildcards, attempt: 57344 + 57344 * attempt,
#         runtime_hrs = lambda wildcards, attempt: 48 * attempt
#     params:
#         sample = SAMPLE,
#         sort_threads = 4,
#         align_threads = int(config['num_cpu_high'] - 4)
#     shell:
#         'winnowmap -W {input.rep_kmer} -t {params.align_threads} -ax map-pb -R "@RG\\tID:1\\tSM:{params.sample}" {input.fasta} {input.reads} | '
#         'samtools sort -m 4096M --threads {params.sort_threads} | samtools view -F 4 -b > {output}'


# rule align_ontec_output_reads:
#     input:
#         reads = expand(
#             'output/ont_ec_subsets/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.clip-ec.geq{{size_fraction}}.fa.gz',
#             zip,
#             filename=ONT_ALL_FILES * len(MBG_KMER_SIZE),
#             kmer=MBG_KMER_SIZE * len(ONT_ALL_FILES),
#             window=MBG_WINDOW_SIZE * len(ONT_ALL_FILES),
#             minscore=[GA_MIN_SCORE] * len(MBG_KMER_SIZE) * len(ONT_ALL_FILES)
#         ),
#         fasta = os.path.join(REFERENCE_FOLDER, '{reference}.fasta'),
#         rep_kmer = os.path.join(REFERENCE_FOLDER, '{reference}.k15.rep-grt09998.txt')
#     output:
#         bam = 'output/read_align/{}_ONTEC_geq{{size_fraction}}_MAP-TO_{{reference}}.k15.bam'.format(SAMPLE)
#     benchmark:
#         'rsrc/output/read_align/{}_ONTEC_geq{{size_fraction}}_MAP-TO_{{reference}}.k15.wmap.rsrc'.format(SAMPLE)
#     log:
#         'log/output/read_align/{}_ONTEC_geq{{size_fraction}}_MAP-TO_{{reference}}.k15.wmap.log'.format(SAMPLE)
#     conda:
#         '../../environment/conda/conda_biotools.yml'
#     threads: config['num_cpu_high']
#     resources:
#         mem_total_mb = lambda wildcards, attempt: 57344 + 57344 * attempt,
#         runtime_hrs = lambda wildcards, attempt: 48 * attempt
#     params:
#         sample = SAMPLE,
#         sort_threads = 4,
#         align_threads = int(config['num_cpu_high'] - 4)
#     shell:
#         'winnowmap -W {input.rep_kmer} -k 15 -t {params.align_threads} -ax map-pb -R "@RG\\tID:1\\tSM:{params.sample}" {input.fasta} {input.reads} | '
#         'samtools sort -m 4096M --threads {params.sort_threads} | samtools view -F 4 -b > {output}'


rule compute_bam_index:
    input:
        '{filepath}/{filename}.bam'
    output:
        '{filepath}/{filename}.bam.bai'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'samtools index -@ {threads} {input}'


rule compute_coverage_per_bp:
    """
    Print positions of non-zero coverage in zero-based coordinates
    """
    input:
        read_ref_aln = 'output/read_align/{read_ref_align}.bam',
        aln_idx = 'output/read_align/{read_ref_align}.bam.bai',
    output:
        'output/read_align_cov/{read_ref_align}.cov.bg.gz',
    benchmark:
        'rsrc/output/read_align_cov/{read_ref_align}.cov.rsrc',
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 + 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        compress_threads = int(config['num_cpu_low'] - 1)
    shell:
        'bedtools genomecov -dz -ibam {input.read_ref_aln} | pigz -p {params.compress_threads} > {output}'


rule cache_positional_coverages:
    input:
        reference = os.path.join(REFERENCE_FOLDER, '{reference}.fasta.fai'),
        coverage = 'output/read_align_cov/{reads}_MAP-TO_{reference}.k{kmer}.cov.bg.gz',
    output:
        cache_ok = 'output/read_align_cov/{reads}_MAP-TO_{reference}.k{kmer}.cache.chk'
    benchmark:
        'rsrc/output/read_align_cov/{reads}_MAP-TO_{reference}.k{kmer}.cache.rsrc'
    resources:
         mem_total_mb = lambda wildcards, attempt: 1024 * attempt
    run:
        import numpy as np
        import gzip as gzip
        chrom_sizes = dict()
        with open(input.reference, 'r') as index:
            for line in index:
                chrom, size = line.split()[:2]
                chrom_sizes[chrom.strip()] = int(size)
        
        current_chrom = None
        chrom_coverage = None
        with gzip.open(input.coverage, 'rt') as track:
            for line in track:
                chrom, position, cov = line.split()
                if chrom != current_chrom:
                    if current_chrom is not None:
                        chrom_cache = output.cache_ok.rsplit('.', 1)[0] + '.{}.npy'.format(current_chrom)
                        np.save(chrom_cache, chrom_coverage, allow_pickle=False)
                    current_chrom = chrom
                    chrom_coverage = np.zeros(chrom_sizes[chrom], dtype=np.int16)
                chrom_coverage[int(position)] = int(cov)
        
        with open(output.cache_ok, 'w') as checkfile:
            pass


rule compute_binned_coverage:
    input:
        reference = os.path.join(REFERENCE_FOLDER, '{reference}.fasta.fai'),
        cache_ok = 'output/read_align_cov/{reads}_MAP-TO_{reference}.k{kmer}.cache.chk'
    output:
        'output/binned_coverage/{reads}_MAP-TO_{reference}.k{kmer}.{binsize}.{chrom}.tsv'
    benchmark:
        'rsrc/output/binned_coverage/{reads}_MAP-TO_{reference}.k{kmer}.{binsize}.{chrom}.poscov.rsrc'
    resources:
         mem_total_mb = lambda wildcards, attempt: 1024 * attempt
    wildcard_constraints:
        chrom='(' + '|'.join(MAIN_CHROMOSOMES) + ')'
    run:
        import numpy as np
        import io as io
        chrom_size = 0
        binsize = int(wildcards.binsize)
        with open(input.reference, 'r') as index:
            for line in index:
                chrom, size = line.split()[:2]
                if chrom == wildcards.chrom:
                    chrom_size = int(size)
                    break
        cache_file = input.cache_ok.rsplit('.', 1)[0] + '.{}.npy'.format(wildcards.chrom)
        try:
            chrom_cov = np.load(cache_file)
        except IOError:
            chrom_cov = np.zeros(chrom_size, dtype=np.int16)
        
        tail_cut = chrom_size // binsize * binsize

        chrom_tail = chrom_cov[tail_cut:]
        chrom_tail.sort()

        chrom_cov = np.reshape(chrom_cov[:tail_cut], (-1, binsize))
        chrom_cov.sort(axis=1)

        cov_bp = (chrom_cov > 0).sum(axis=1)
        cov_pct = (cov_bp / binsize * 100).round(2)
        summed_cov = chrom_cov.sum(axis=1)
        mean_cov = chrom_cov.mean(axis=1).round(2)
        median_cov = chrom_cov[:, binsize // 2]
        mean_nzcov = (summed_cov / cov_bp).round(2)
        np.nan_to_num(mean_nzcov, copy=False, nan=0.0, posinf=1000000000, neginf=-1)

        out_buffer = io.StringIO()
        _ = out_buffer.write(
            '\t'.join([
                '#chrom',
                'start',
                'end',
                'mean_cov',
                'median_cov',
                'mean_nzcov',
                'cov_bp',
                'cov_bp_pct',
                'summed_cov'
            ]) + '\n'
        )
        chrom = wildcards.chrom
        row_template = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
        for i, start in enumerate(range(0, tail_cut, binsize)):
            end = start + binsize
            _ = out_buffer.write(
                row_template.format(
                    chrom,
                    start,
                    end,
                    mean_cov[i],
                    median_cov[i],
                    mean_nzcov[i],
                    cov_bp[i],
                    cov_pct[i],
                    summed_cov[i]
                )
            )
        
        tail_mean_nzcov = (chrom_tail.sum() / (chrom_tail > 0).sum()).round(2)
        if np.isnan(tail_mean_nzcov):
            tail_mean_nzcov = 0.

        _ = out_buffer.write(
                row_template.format(
                    chrom,
                    tail_cut,
                    tail_cut + chrom_tail.size,
                    chrom_tail.mean().round(2),
                    chrom_tail[chrom_tail.size // 2],
                    tail_mean_nzcov,
                    (chrom_tail > 0).sum(),
                    ((chrom_tail > 0).sum() / chrom_tail.size * 100).round(2),
                    chrom_tail.sum()
                )
            )

        with open(output[0], 'w') as table:
            _ = table.write(out_buffer.getvalue())
    # END OF RUN BLOCK


rule merge_binned_coverage:
    input:
        tables = expand('output/binned_coverage/{{reads}}_MAP-TO_{{reference}}.k{{kmer}}.{{binsize}}.{chrom}.tsv',
                        chrom=MAIN_CHROMOSOMES)
    output:
        table = 'output/binned_coverage/{reads}_MAP-TO_{reference}.k{kmer}.{binsize}.merged.tsv'
    run:
        import pandas as pd
        chrom_coverages = []
        for tsv_file in input.tables:
            df = pd.read_csv(tsv_file, sep='\t', header=0, index_col=None)
            chrom_coverages.append(df)
        chrom_coverages = pd.concat(
            chrom_coverages,
            axis=0,
            ignore_index=False
        )
        chrom_coverages.sort_values(['#chrom', 'start', 'end'], ascending=True, inplace=True)
        chrom_coverages.to_csv(
            output.table,
            sep='\t',
            header=True,
            index=False
        )


# rule merge_short_reads:
#     input:
#         reads1 = '{}_1_val_1.cor.fq.gz'.format(SHORT_READS_PATH),
#         reads2 = '{}_2_val_2.cor.fq.gz'.format(SHORT_READS_PATH)
#     output:
#         fastq = 'output/fastq/{sample}_SHORT.fastq.gz'
#     benchmark:
#         'rsrc/output/fastq/{sample}_SHORT.mrg.rsrc'
#     conda:
#         '../../environment/conda/conda_biotools.yml'
#     threads: config['num_cpu_medium']
#     resources:
#         runtime_hrs = lambda wildcards, attempt: 6 * attempt
#     params:
#         threads = int(config['num_cpu_medium'] // 2)
#     shell:
#         'pigz -p {params.threads} -d -c {input.reads1} {input.reads2} | '
#         'pigz -p {params.threads} --best > {output.fastq}'


# rule merge_ontec_reads:
#     """
#     This rule merges ONTEC reads per MBG and minscore parameterization
#     Goal: one file of input reads per combination of k/w/minscore for Bifrost
#     """
#     input:
#         reads = expand(
#             'output/ont_ec_subsets/{filename}_MAP-TO_mbg-k{{kmer}}-w{{window}}.ms{{minscore}}.clip-ec.geq0.fa.gz',
#             filename=ONT_ALL_FILES * len(MBG_KMER_SIZE),
#         )
#     output:
#         fastq = 'output/fastq/{}_ONTEC_mbg-k{{kmer}}-w{{window}}.ms{{minscore}}.geq0.fastq.gz'
#     benchmark:
#         'rsrc/output/fastq/{}_ONTEC_mbg-k{{kmer}}-w{{window}}.ms{{minscore}}.geq0.mrg.rsrc'
#     conda:
#         '../../environment/conda/conda_biotools.yml'
#     threads: config['num_cpu_medium']
#     resources:
#         runtime_hrs = lambda wildcards, attempt: 6 * attempt
#     params:
#         threads = int(config['num_cpu_medium'] // 2)
#     shell:
#         'pigz -p {params.threads} -d -c {input.reads} | '
#         'pigz -p {params.threads} --best > {output.fastq}'


# rule hifiasm_ontec_only_assembly:
#     """
#     """
#     input:
#         container = 'hifiasm-v0142r315.sif',
#         ontec_reads = expand(
#             'output/ont_ec_subsets/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.clip-ec.geq{{size_fraction}}.fa.gz',
#             zip,
#             filename=ONT_ALL_FILES * len(MBG_KMER_SIZE),
#             kmer=MBG_KMER_SIZE * len(ONT_ALL_FILES),
#             window=MBG_WINDOW_SIZE * len(ONT_ALL_FILES),
#             minscore=[GA_MIN_SCORE] * len(MBG_KMER_SIZE) * len(ONT_ALL_FILES)
#         )
#     output:
#         primary_unitigs = 'output/assembly/layout/{}_ontec-only_geq{{size_fraction}}/{}_ontec-only_geq{{size_fraction}}.p_utg.gfa'.format(SAMPLE, SAMPLE),
#         primary_contigs = 'output/assembly/layout/{}_ontec-only_geq{{size_fraction}}/{}_ontec-only_geq{{size_fraction}}.p_ctg.gfa'.format(SAMPLE, SAMPLE),
#         raw_unitigs = 'output/assembly/layout/{}_ontec-only_geq{{size_fraction}}/{}_ontec-only_geq{{size_fraction}}.r_utg.gfa'.format(SAMPLE, SAMPLE),
#         discard = multiext(
#             'output/assembly/layout/{}_ontec-only_geq{{size_fraction}}/{}_ontec-only_geq{{size_fraction}}'.format(SAMPLE, SAMPLE),
#             '.a_ctg.gfa', '.a_ctg.noseq.gfa',
#             '.ec.bin', '.ovlp.reverse.bin', '.ovlp.source.bin',
#             '.p_ctg.noseq.gfa', '.p_utg.noseq.gfa',
#             '.r_utg.noseq.gfa'
#         )
#     log:
#         hifiasm = 'log/output/assembly/layout/{}_ontec-only_{{size_fraction}}.hifiasm.log'.format(SAMPLE)
#     benchmark:
#         'rsrc/output/assembly/layout/{}_ontec-only_{{size_fraction}}.hifiasm'.format(SAMPLE) + '.t{}.rsrc'.format(config['num_cpu_max']) 
# #    conda:
# #        '../environment/conda/conda_biotools.yml'       
#     threads: config['num_cpu_max']
#     resources:
#         mem_total_mb = lambda wildcards, attempt: set_hifiasm_memory(wildcards, attempt, RUN_SYSTEM),
#         runtime_hrs = lambda wildcards, attempt: 52 + 52 * attempt
#     params:
#         prefix = lambda wildcards, output: output.primary_contigs.rsplit('.', 2)[0],
#         singularity = '' if not config.get('env_module_singularity', False) else 'module load {} ; '.format(config['env_module_singularity'])
#     shell:
#         '{params.singularity} singularity exec {input.container} '
#         'hifiasm -o {params.prefix} -t {threads} {input.ontec_reads} &> {log.hifiasm}'


# rule write_reads_fofn:
#     input:
#         short_reads = 'output/fastq/{sample}_SHORT.fastq.gz',
#         ontec_reads = expand(
#             'output/fastq/{{sample}}_ONTEC_mbg-k{kmer}-w{window}.ms{minscore}.geq0.fastq.gz',
#             zip,
#             kmer=MBG_KMER_SIZE,
#             window=MBG_WINDOW_SIZE,
#             minscore=[GA_MIN_SCORE] * len(MBG_KMER_SIZE)
#         ),
#         hifi_reads = HIFI_READS_PATH
#     output:
#         fofn = 'output/fofn/{sample}_all_reads.fofn'
#     run:
#         import os
#         with open(output.fofn, 'w') as listing:
#             for entry in input:
#                 assert os.path.isfile(entry), 'Error record: {} / {}'.format(entry, input)
#                 full_path = os.path.abspath(entry)
#                 _ = listing.write(full_path + '\n')


# if RUN_SYSTEM == 'valet':
#     bifrost_cpu = config['num_cpu_max']
# else:
#     bifrost_cpu = config['num_cpu_high']


# rule build_colored_dbg:
#     input:
#         fofn = 'output/fofn/{sample}_all_reads.fofn'
#     output:
#         gfa = 'output/cdbg/{sample}.{strategy}.gfa',
#         colors = 'output/cdbg/{sample}.{strategy}.bfg_colors'
#     log:
#        'log/output/cdbg/{sample}.{strategy}.build.log',
#     benchmark:
#         'rsrc/output/cdbg/{sample}.{strategy}.build.rsrc',
#     wildcard_constraints:
#         strategy = '(all\-ref|all\-seq)'
#     conda: '../../environment/conda/conda_biotools.yml'
#     threads: bifrost_cpu
#     resources:
#         mem_total_mb = lambda wildcards, attempt: 63488 * attempt,
#         runtime_hrs = lambda wildcards, attempt: 48 * attempt
#     params:
#         kmer_size = 31,
#         out_prefix = lambda wildcards, output: output.gfa.rsplit('.', 1)[0],
#         input_type = lambda wildcards: '--input-ref-file' if wildcards.strategy == 'all-ref' else '--input-seq-file'
#     shell:
#         'Bifrost build {params.input_type} {input.fofn} '
#         '--output-file {params.out_prefix} --threads {threads} --colors --kmer-length {params.kmer_size} '
#         '--verbose &> {log}'


# rule count_cdbg_kmers:
#     input:
#         container = 'cdbg_kmc.sif',
#         gfa = 'output/cdbg/{sample}.{strategy}.gfa',
#         colors = 'output/cdbg/{sample}.{strategy}.bfg_colors'
#     output:
#         'output/cdbg_kcount/{sample}.{strategy}.kmc.txt'
#     log:
#        'log/output/cdbg_kcount/{sample}.{strategy}.kmc.log',
#     benchmark:
#         'rsrc/output/cdbg_kcount/{sample}.{strategy}.kmc.rsrc',
#     threads: config['num_cpu_low']
#     resources:
#         mem_total_mb = lambda wildcards, attempt: 172032 * attempt,
#         runtime_hrs = lambda wildcards, attempt: 36 * attempt
#     params:
#         kmer_size = 31,
#         singularity = '' if not config.get('env_module_singularity', False) else 'module load {} ; '.format(config['env_module_singularity'])
#     shell:
#         '{params.singularity} singularity exec {input.container} venn_diagram '
#             '{input.gfa} {input.colors} {params.kmer_size} {threads} {output} &> {log}'


# # rule dump_mbg_tigs_to_fasta:
# #     input:
# #         gfa = 'output/mbg_hifi_clean/{sample}_HiFi.mbg-k{kmer}-w{window}.clean.gfa'
# #     output:
# #         fasta = 'output/graph_sequences/{sample}_HiFi_MBG.k{kmer}-w{window}.fasta',
# #         stats = 'output/stats/{sample}_HiFi_MBG.k{kmer}-w{window}.contig.stats',
# #         mapping = 'output/stats/{sample}_HiFi_MBG.k{kmer}-w{window}.read-contig.map'
# #     log:
# #         'log/output/sequences/{sample}_HiFi_MBG.k{kmer}-w{window}.dump.log'
# #     benchmark:
# #         'rsrc/output/sequences/{sample}_HiFi_MBG.k{kmer}-w{window}.dump.rsrc'
# #     conda:
# #         '../../environment/conda/conda_pyscript.yml'
# #     threads: config['num_cpu_low']
# #     resources:
# #         mem_per_cpu_mb = lambda wildcards, attempt: 24768 * attempt,
# #         mem_total_mb = lambda wildcards, attempt: 24768 * attempt,
# #         runtime_hrs = lambda wildcards, attempt: attempt * attempt
# #     params:
# #         script_exec = lambda wildcards: find_script_path('gfa_to_fasta.py')
# #     shell:
# #         '{params.script_exec} --gfa {input.gfa} --n-cpus {threads} '
# #         '--out-fasta {output.fasta} --out-map {output.mapping} '
# #         '--out-stats {output.stats} &> {log}'


# rule dump_unitgs_to_fasta:
#     input:
#         gfa = 'output/assembly/layout/{assembly_graph}/{assembly_graph}.{tigs}.gfa'
#     output:
#         fasta = 'output/graph_sequences/{assembly_graph}.{tigs}.fasta',
#         stats = 'output/stats/{assembly_graph}.{tigs}.contig.stats',
#         mapping = 'output/stats/{assembly_graph}.{tigs}.read-contig.map'
#     log:
#         'log/output/sequences/{assembly_graph}.{tigs}.dump.log'
#     benchmark:
#         'rsrc/output/sequences/{assembly_graph}.{tigs}.dump.rsrc'
#     conda:
#         '../../environment/conda/conda_pyscript.yml'
#     wildcard_constraints:
#         tigs = '(r_utg|p_ctg|p_utg)'
#     threads: config['num_cpu_low']
#     resources:
#         mem_per_cpu_mb = lambda wildcards, attempt: 16384 * attempt,
#         mem_total_mb = lambda wildcards, attempt: 16384 * attempt,
#         runtime_hrs = lambda wildcards, attempt: attempt * attempt
#     params:
#         script_exec = lambda wildcards: find_script_path('gfa_to_fasta.py')
#     shell:
#         '{params.script_exec} --gfa {input.gfa} --n-cpus {threads} '
#         '--out-fasta {output.fasta} --out-map {output.mapping} '
#         '--out-stats {output.stats} &> {log}'


# # rule align_tig_sequences_to_reference:
# #     """
# #     winnowmap tends to segfaults for unknown reasons
# #     on certain tig datasets - switch to minimap
# #     """
# #     input:
# #         reference = os.path.join(REFERENCE_FOLDER, '{reference}.fasta'),
# #         rep_kmer = os.path.join(REFERENCE_FOLDER, '{reference}.k19.rep-grt09998.txt'),
# #         tig_seq = 'output/graph_sequences/{assembly_graph}.{tigs}.fasta',
# #     output:
# #         bam = 'output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k19.wmap.bam'
# #     benchmark:
# #         'rsrc/output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k19.wmap.rsrc'
# #     log:
# #         'log/output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k19.wmap.log'
# #     conda:
# #         '../../environment/conda/conda_biotools.yml'
# #     threads: config['num_cpu_high']
# #     resources:
# #         mem_total_mb = lambda wildcards, attempt: 36864 + 36864 * attempt,
# #         runtime_hrs = lambda wildcards, attempt: 23 * attempt
# #     params:
# #         sample = lambda wildcards: wildcards.assembly_graph.split('_')[0],
# #         sort_threads = 4,
# #         align_threads = int(config['num_cpu_high'] - 4)
# #     shell:
# #         'winnowmap -W {input.rep_kmer} -k 19 -t {params.align_threads} -ax map-pb -R "@RG\\tID:1\\tSM:{params.sample}" {input.reference} {input.tig_seq} | '
# #         'samtools sort -m 4096M --threads {params.sort_threads} --output-fmt BAM -l 9 > {output.bam}'


# rule align_tig_sequences_to_reference:
#     """
#     winnowmap tends to segfaults for unknown reasons
#     on certain tig datasets - switch to minimap
#     """
#     input:
#         reference = os.path.join(REFERENCE_FOLDER, '{reference}.fasta'),
#         tig_seq = 'output/graph_sequences/{assembly_graph}.{tigs}.fasta',
#     output:
#         bam = 'output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k19.mmap.bam'
#     benchmark:
#         'rsrc/output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k19.mmap.rsrc'
#     log:
#         'log/output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k19.mmap.log'
#     conda:
#         '../../environment/conda/conda_biotools.yml'
#     threads: config['num_cpu_high']
#     resources:
#         mem_total_mb = lambda wildcards, attempt: 36864 + 36864 * attempt,
#         runtime_hrs = lambda wildcards, attempt: 23 * attempt
#     params:
#         sample = lambda wildcards: wildcards.assembly_graph.split('_')[0],
#         sort_threads = 4,
#         align_threads = int(config['num_cpu_high'] - 4)
#     shell:
#         'minimap2 --secondary=no --eqx -Y -ax asm5 '
#         '-t {threads} -R "@RG\\tID:1\\tSM:{params.sample}" {input.reference} {input.tig_seq} | '
#         'samtools sort -m 4096M --threads {params.sort_threads} --output-fmt BAM -l 9 > {output.bam}'


# rule dump_tig_alignments_to_bed:
#     input:
#         bam = 'output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k{kmer}.mmap.bam'
#     output:
#         bed = 'output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k{kmer}.bed.gz'
#     conda:
#         '../../environment/conda/conda_biotools.yml'
#     threads: config['num_cpu_low']
#     resources:
#         mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
#         runtime_hrs = lambda wildcards, attempt: 6 * attempt
#     params:
#         compress_threads = int(config['num_cpu_low'] - 1)
#     shell:
#         'bedtools bamtobed -splitD -i {input.bam} | pigz -p {params.compress_threads} > {output.bed}'


# rule tig_alignments_in_roi:
#     input:
#         bed_aln = 'output/alignments/tig_to_ref/{assembly_graph}.{tigs}_MAP-TO_{reference}.k{kmer}.bed.gz',
#         bed_regions = os.path.join(REFERENCE_FOLDER, 'T2Tv1_38p13Y_chm13.{subset}.bed')
#     output:
#         'output/tig_intersect/tig_in_roi/{assembly_graph}.{tigs}_MAP-TO_{reference}.k{kmer}.{subset}.tsv'
#     conda:
#         '../../environment/conda/conda_biotools.yml'
#     shell:
#         'bedtools intersect -wo -a {input.bed_aln} -b {input.bed_regions} > {output}'


# rule subset_assembly_graph:
#     input:
#         table = 'output/tig_intersect/tig_in_roi/{assembly_graph}.{tigs}_MAP-TO_{reference}.k19.{subset}.tsv',
#         graph = 'output/assembly/layout/{assembly_graph}/{assembly_graph}.{tigs}.gfa'
#     output:
#         graph = 'output/subset_graphs/{assembly_graph}.{tigs}_MAP-TO_{reference}.k19.{subset}.gfa',
#         table = 'output/subset_graphs/{assembly_graph}.{tigs}_MAP-TO_{reference}.k19.{subset}.csv',
#     log:
#         'log/output/subset_graphs/{assembly_graph}.{tigs}_MAP-TO_{reference}.k19.{subset}.subset.log',
#     conda:
#         '../../environment/conda/conda_pyscript.yml'
#     params:
#         script_exec = lambda wildcards: find_script_path('gfa_subset.py')
#     shell:
#         '{params.script_exec} --debug --input-gfa {input.graph} --input-table {input.table} '
#         '--output-gfa {output.graph} --output-table {output.table} &> {log}'
