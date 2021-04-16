localrules: master

# config values data
SAMPLE = config['sample']
RUN_SYSTEM = config['system']

ONT_HPG_FILES = config['ont_hpg_files'][SAMPLE]
ONT_ALL_FILES = ONT_HPG_FILES

ONT_RAW_FOLDER = config['ont_raw_folder'][SAMPLE][RUN_SYSTEM]
HIFI_READS_PATH = config['hifi_reads_path'][SAMPLE][RUN_SYSTEM]

# config values system
RUNTIME_GA = config['runtime']['GraphAligner'][RUN_SYSTEM]
RUNTIME_ES = config['runtime']['ExtractSequence'][RUN_SYSTEM]

MEMORY_MBG_LOW = config['memory']['MBG_low'][SAMPLE]  # large k-mer
MEMORY_MBG_HIGH = config['memory']['MBG_high'][SAMPLE]  # small k-mer

MEMORY_GA = config['memory']['GraphAligner'][SAMPLE]

MEMORY_HIFIASM = config['memory']['HifiAsm']

# config values parameters
MBG_KMER_SIZE = config['params']['MBG']['kmer']
MBG_WINDOW_SIZE = config['params']['MBG']['window']

GA_MIN_SCORE = config['params']['GraphAligner']['minscore']

# size fraction zero = dump all reads
READ_SIZE_FRACTIONS = config['params']['SizeFractions']


if SAMPLE == 'HG00733':

    # ===== DEBUG
    MBG_KMER_SIZE = [k for k in MBG_KMER_SIZE if int(k) > 1000]
    MBG_WINDOW_SIZE = [w for w in MBG_WINDOW_SIZE if int(w) > 500]
    # ===== DEBUG

    ONT_SPLIT_FILES = config['ont_split_files']
    ONT_ALL_FILES.append('HG00733_EEE_ONT')
    rule merge_ont_splits:
        input:
            fastq = expand('{}/{{filename}}.fastq.gz'.format(ONT_RAW_FOLDER), filename=ONT_SPLIT_FILES)
        output:
            fasta = 'input/ont/HG00733_EEE_ONT.fa.gz'
        benchmark:
            'rsrc/input/ont/HG00733_EEE_ONT.merge.rsrc'
        conda:
            '../../environment/conda/conda_biotools.yml'
        threads: config['num_cpu_low']
        resources:
            mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
            runtime_hrs = lambda wildcards, attempt: attempt ** attempt
        params:
            pigz_cpu = lambda wildcards: int(config['num_cpu_low']) // 2
        shell:
            'pigz -p {params.pigz_cpu} -c -d {input.fastq} | '
            'seqtk seq -S -A -l 120 -C | '
            'pigz -p {params.pigz_cpu} > {output.fasta}'


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


rule clean_hpg_ont:
    input:
        fastq = '{}/{{filename}}.fastq.gz'.format(ONT_RAW_FOLDER)
    output:
        fasta = 'input/ont/{filename}.fa.gz'
    benchmark:
        'rsrc/input/ont/{filename}.clean.rsrc'
    wildcard_constraints:
        filename = '(' + '|'.join(ONT_HPG_FILES) + ')'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt ** attempt
    params:
        pigz_cpu = lambda wildcards: int(config['num_cpu_low']) // 2
    shell:
        'pigz -p {params.pigz_cpu} -c -d {input.fastq} | '
        'seqtk seq -S -A -l 120 -C {input.fastq} | '
        'pigz -p {params.pigz_cpu} > {output.fasta}'


def set_mbg_memory(wildcards, attempt):

    if int(wildcards.kmer) < 1000:
        return MEMORY_MBG_HIGH + MEMORY_MBG_HIGH * attempt
    else:
        return MEMORY_MBG_LOW + MEMORY_MBG_LOW * attempt


rule build_hifi_read_graph:
    input:
        HIFI_READS_PATH
    output:
        'output/mbg_hifi/{}_HiFi.mbg-k{{kmer}}-w{{window}}.gfa'.format(SAMPLE)
    log:
        'log/output/mbg_hifi/{}_HiFi.mbg-k{{kmer}}-w{{window}}.log'.format(SAMPLE)
    benchmark:
        'rsrc/output/mbg_hifi/{}_HiFi.mbg-k{{kmer}}-w{{window}}.rsrc'.format(SAMPLE)
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = set_mbg_memory,
        runtime_hrs = lambda wildcards, attempt: attempt + attempt * attempt
    shell:
        'MBG -i {input} -o {output} -t {threads} --blunt -k {wildcards.kmer} -w {wildcards.window} &> {log}'


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


rule strip_sequences_from_graph:
    input:
        gfa = 'output/mbg_hifi_clean/{sample}_HiFi.mbg-k{kmer}-w{window}.clean.gfa'
    output:
        gfa = 'output/mbg_hifi_clean/{sample}_HiFi.mbg-k{kmer}-w{window}.clean.noseq.gfa'
    benchmark:
        'rsrc/output/mbg_hifi_clean/{sample}_HiFi.mbg-k{kmer}-w{window}.noseq.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'gfatools view -S {input.gfa} > {output.gfa}'


rule ont_error_correction:
    """
    """
    input:
        graph = 'output/mbg_hifi_clean/{}_HiFi.mbg-k{{kmer}}-w{{window}}.clean.gfa'.format(SAMPLE),
        reads = 'input/ont/{filename}.fa.gz',
    output:
        gaf = 'output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.gaf',
        ec_reads = 'output/ont_ec/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.clip-ec.fa.gz',
    log:
        'log/output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.ga.log'
    benchmark:
        'rsrc/output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.ga.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: MEMORY_GA * attempt,
        runtime_hrs = lambda wildcards, attempt: RUNTIME_GA * attempt * attempt
    shell:
        'GraphAligner -t {threads} -g {input.graph} -f {input.reads} '
            '-x dbg --min-alignment-score {wildcards.minscore} '
            '--corrected-clipped-out {output.ec_reads} '
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
        tsv = expand('output/seq_stats/input/ont/{filename}.stats.tsv.gz', filename=ONT_ALL_FILES)
    output:
        summary = 'output/seq_stats/summary/{}_ONT_input.stats.tsv'.format(SAMPLE),
        dump = 'output/seq_stats/summary/{}_ONT_input.stats.pck'.format(SAMPLE)
    log:
        'log/output/seq_stats/summary/{}_ONT_input.stats.log'.format(SAMPLE)
    benchmark:
        'rsrc/output/seq_stats/summary/{}_ONT_input.stats.rsrc'.format(SAMPLE)
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


rule extract_ec_reads_by_size:
    input:
        ec_reads = 'output/ont_ec/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.clip-ec.fa.gz',
        read_list = 'output/read_info/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.geq{size_fraction}.read-ec.txt'
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
        'seqtk subseq -l 120 /dev/stdin {input.read_list} | '
        'pigz -p {params.pigz_cpu} --best > {output}'


rule compute_stats_corrected_reads:
    input:
        tsv = expand(
            'output/seq_stats/output/ont_ec_subsets/{filename}_MAP-TO_mbg-k{{kmer}}-w{{window}}.ms{{minscore}}.clip-ec.geq0.stats.tsv.gz',
            filename=ONT_ALL_FILES
            )
    output:
        summary = 'output/seq_stats/summary/{}_ONT_clip-ec.k{{kmer}}-w{{window}}.ms{{minscore}}.stats.tsv'.format(SAMPLE),
        dump = 'output/seq_stats/summary/{}_ONT_clip-ec.k{{kmer}}-w{{window}}.ms{{minscore}}.stats.pck'.format(SAMPLE)
    log:
        'log/output/seq_stats/summary/{}_ONT_clip-ec.k{{kmer}}-w{{window}}.ms{{minscore}}.stats.log'.format(SAMPLE)
    benchmark:
        'rsrc/output/seq_stats/summary/{}_ONT_clip-ec.k{{kmer}}-w{{window}}.ms{{minscore}}.stats.rsrc'.format(SAMPLE)
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

    if run_system == 'valet':
        return MEMORY_HIFIASM
    else:
        return MEMORY_HIFIASM * attempt


rule hifiasm_hifi_ontec_assembly:
    """
    """
    input:
        container = 'hifiasm-v0142r315.sif',
        hifi_reads = os.path.basename(HIFI_READS_PATH),
        ontec_reads = expand(
            'output/ont_ec_subsets/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.clip-ec.geq{{size_fraction}}.fa.gz',
            zip,
            filename=ONT_ALL_FILES * len(MBG_KMER_SIZE),
            kmer=MBG_KMER_SIZE * len(ONT_ALL_FILES),
            window=MBG_WINDOW_SIZE * len(ONT_ALL_FILES),
            minscore=[GA_MIN_SCORE] * len(MBG_KMER_SIZE) * len(ONT_ALL_FILES)
        )
    output:
        primary_unitigs = 'output/assembly/layout/{}_hifi-ontec_{{size_fraction}}/{}_hifi-ontec_geq{{size_fraction}}.p_utg.gfa'.format(SAMPLE, SAMPLE),
        primary_contigs = 'output/assembly/layout/{}_hifi-ontec_{{size_fraction}}/{}_hifi-ontec_geq{{size_fraction}}.p_ctg.gfa'.format(SAMPLE, SAMPLE),
        raw_unitigs = 'output/assembly/layout/{}_hifi-ontec_{{size_fraction}}/{}_hifi-ontec_geq{{size_fraction}}.r_utg.gfa'.format(SAMPLE, SAMPLE),
        discard = multiext(
            'output/assembly/layout/{}_hifi-ontec_{{size_fraction}}/{}_hifi-ontec_geq{{size_fraction}}'.format(SAMPLE, SAMPLE),
            '.a_ctg.gfa', '.a_ctg.noseq.gfa',
            '.ec.bin', '.ovlp.reverse.bin', '.ovlp.source.bin',
            '.p_ctg.noseq.gfa', '.p_utg.noseq.gfa',
            '.r_utg.noseq.gfa'
        )
    log:
        hifiasm = 'log/output/assembly/layout/{}_hifi-ontec_{{size_fraction}}.hifiasm.log'.format(SAMPLE)
    benchmark:
        'rsrc/output/assembly/layout/{}_hifi-ontec_{{size_fraction}}.hifiasm'.format(SAMPLE) + '.t{}.rsrc'.format(config['num_cpu_max']) 
#    conda:
#        '../environment/conda/conda_biotools.yml'       
    threads: config['num_cpu_max']
    resources:
        mem_total_mb = lambda wildcards, attempt: set_hifiasm_memory(wildcards, attempt, RUN_SYSTEM),
        runtime_hrs = lambda wildcards, attempt: 48 + 48 * attempt
    params:
        prefix = lambda wildcards, output: output.primary_contigs.rsplit('.', 2)[0],
        singularity = '' if not config.get('env_module_singularity', False) else 'module load {} ; '.format(config['env_module_singularity'])
    shell:
        '{params.singularity} singularity exec {input.container} '
        'hifiasm -o {params.prefix} -t {threads} {input.hifi_reads} {input.ontec_reads} &> {log.hifiasm}'


rule master:
    input:
        expand(
            rules.strip_sequences_from_graph.output.gfa,
            zip,
            sample=[SAMPLE] * len(MBG_KMER_SIZE),
            kmer=MBG_KMER_SIZE,
            window=MBG_WINDOW_SIZE
        ),
        rules.compute_stats_input_reads.output.summary,
        expand(
            rules.compute_stats_corrected_reads.output.summary,
            zip,
            kmer=MBG_KMER_SIZE,
            window=MBG_WINDOW_SIZE,
            minscore=[GA_MIN_SCORE] * len(MBG_KMER_SIZE)
        ),
        expand(
            rules.hifiasm_hifi_ontec_assembly.output.raw_unitigs,
            size_fraction=[0, 100000]
        )
