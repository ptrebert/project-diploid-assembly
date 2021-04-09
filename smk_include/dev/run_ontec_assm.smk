localrules: master

### INPUTFILES

eee_ont = [
    'HG00733_ont-ul_SRR13356541_1.fastq.gz',
    'HG00733_ont-ul_SRR13362355_2.fastq.gz',
    'HG00733_ont-ul_SRR13362356_3.fastq.gz',
    'HG00733_ont-ul_SRR13362357_4.fastq.gz',
    'HG00733_ont-ul_SRR13362358_5.fastq.gz',
    'HG00733_ont-ul_SRR13362359_6.fastq.gz'
]

hpg_ont = [
    'HG00733_1_Guppy_4.2.2_prom.fastq.gz',
    'HG00733_2_Guppy_4.2.2_prom.fastq.gz',
    'HG00733_3_Guppy_4.2.2_prom.fastq.gz'
]

### OUTPUTFILES

ont_read_files = [
    'HG00733_EEE_ONT',
    'HG00733_1_Guppy_4.2.2_prom',
    'HG00733_2_Guppy_4.2.2_prom',
    'HG00733_3_Guppy_4.2.2_prom'
]

output_files = [
    'input/ont/HG00733_EEE_ONT.fa.gz',
    'output/mbg_hifi/HG00733_HiFi.mbg-k5001-w2000.noseq.gfa',
    'output/mbg_hifi/HG00733_HiFi.mbg-k2001-w1000.noseq.gfa',
    'output/mbg_hifi/HG00733_HiFi.mbg-k501-w100.noseq.gfa',
]

pattern_reads = 'output/ont_ec/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ms{minscore}.clip-ec.fa.gz'
pattern_stats = 'output/seq_stats/summary/HG00733_ONT_clip-ec.k{kmer}-w{window}.ms{minscore}.stats.tsv'

for orf in ont_read_files:
    for k, w in [(5001,2000), (2001,1000), (501,100)]:        
        for ms in [0]:
            tmp = pattern_reads.format(**values)
            output_files.append(tmp)
            values = {
                'filename': orf,
                'kmer': k,
                'window': w,
                'minscore': ms
            }
            tmp = pattern_stats.format(**values)
            output_files.append(tmp)


rule merge_eee_ont:
    input:
        fastq = expand('/gpfs/project/projects/medbioinf/data/hg00733_ont/{filename}', filename=eee_ont)
    output:
        fasta = 'input/ont/HG00733_EEE_ONT.fa.gz'
    benchmark:
        'rsrc/input/ont/HG00733_EEE_ONT.merge.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(2048 * attempt // config['num_cpu_low']),
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt ** attempt
    params:
        pigz_cpu = lambda wildcards: int(config['num_cpu_low']) // 2
    shell:
        'pigz -p {params.pigz_cpu} -c -d {input.fastq} | seqtk seq -S -A -l 120 -C | pigz -p {params.pigz_cpu} > {output.fasta}'


rule clean_hpg_ont:
    input:
        fastq = '/gpfs/project/projects/medbioinf/data/hg00733_ont/{filename}.fastq.gz'
    output:
        fasta = 'input/ont/{filename}.fa.gz'
    benchmark:
        'rsrc/input/ont/{filename}.clean.rsrc'
    wildcard_constraints:
        filename = '(' + '|'.join([f.rsplit('.', 2)[0] for f in hpg_ont]) + ')'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(2048 * attempt // config['num_cpu_low']),
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt + attempt ** attempt
    params:
        pigz_cpu = lambda wildcards: int(config['num_cpu_low']) - 1
    shell:
        'seqtk seq -S -A -l 120 -C {input.fastq} | pigz -p {params.pigz_cpu} > {output.fasta}'


rule build_hifi_read_graph:
    input:
        '/gpfs/project/projects/medbioinf/projects/hifi_v13/run_folder/input/fastq/HG00733_hgsvc_pbsq2-ccs_1000.fastq.gz'
    output:
        'output/mbg_hifi/HG00733_HiFi.mbg-k{kmer}-w{window}.gfa'
    log:
        'log/output/mbg_hifi/HG00733_HiFi.mbg-k{kmer}-w{window}.log'
    benchmark:
        'rsrc/output/mbg_hifi/HG00733_HiFi.mbg-k{kmer}-w{window}.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(245760 + (122880 if int(wildcards.kmer) > 5000 else 245760) * attempt // config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 245760 + (122880 if int(wildcards.kmer) > 5000 else 245760) * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    shell:
        'MBG -i {input} -o {output} -t {threads} --blunt -k {wildcards.kmer} -w {wildcards.window} &> {log}'


rule clean_mbg_graph:
    """
    https://github.com/maickrau/MBG/issues/1

    ### I've used vg: vg view -Fv graph.gfa | vg mod -n -U 100 - | vg view - > blunt-graph.gfa
    """
    input:
        'output/mbg_hifi/HG00733_HiFi.mbg-k{kmer}-w{window}.gfa'
    output:
        'output/mbg_hifi_clean/HG00733_HiFi.mbg-k{kmer}-w{window}.clean.gfa'
    log:
        'log/output/mbg_hifi_clean/HG00733_HiFi.mbg-k{kmer}-w{window}.clean.log'
    benchmark:
        'rsrc/output/mbg_hifi_clean/HG00733_HiFi.mbg-k{kmer}-w{window}.clean.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 16384 * attempt,
        mem_total_mb = lambda wildcards, attempt: 16384 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt ** attempt
    shell:
        'vg view -Fv {input} | vg mod -n -U 100 - 2> {log} | vg view - > {output}'


rule strip_sequences_from_graph:
    input:
        'output/mbg_hifi_clean/HG00733_HiFi.mbg-k{kmer}-w{window}.clean.gfa'
    output:
        'output/mbg_hifi_clean/HG00733_HiFi.mbg-k{kmer}-w{window}.clean.noseq.gfa'
    benchmark:
        'rsrc/output/mbg_hifi_clean/HG00733_HiFi.mbg-k{kmer}-w{window}.noseq.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt ** attempt
    shell:
        'gfatools view -S {input} > {output}'


rule ont_error_correction:
    input:
        graph = 'output/mbg_hifi_clean/HG00733_HiFi.mbg-k{kmer}-w{window}.clean.gfa',
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
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(32768 + 32768 * attempt // config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: 8 * attempt
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
        mem_per_cpu_mb = lambda wildcards, attempt: 4096,
        mem_total_mb = lambda wildcards, attempt: 4096,
        runtime_hrs = lambda wildcards, attempt: 12 * attempt
    params:
        threads = int(config['num_cpu_low'] - 1)
    shell:
        'seqtk comp {input} | pigz -p {params.threads} --best > {output}'


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


rule compute_stats_input_reads:
    input:
        fastq = expand(
            'output/seq_stats/input/ont/{filename}.stats.tsv.gz',
            filename=ont_read_files
        )
    output:
        summary = 'output/seq_stats/summary/HG00733_ONT_input.stats.tsv',
        dump = 'output/seq_stats/summary/HG00733_ONT_input.stats.pck'
    log:
        'log/output/seq_stats/summary/HG00733_ONT_input.stats.log'
    benchmark:
        'rsrc/output/seq_stats/summary/HG00733_ONT_input.stats.rsrc'
    conda:
        '../../environment/conda/conda_pyscript.yml'
    threads: config['num_cpu_low']
    params:
        script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
        '{params.script_exec} --debug --input-files {input.fastq} '
        '--output {output.dump} --summary-output {output.summary} '
        '--num-cpu {threads} --genome-size 3100000000 &> {log}'


rule compute_stats_corrected_reads:
    input:
        fasta = expand(
            'output/seq_stats/output/ont_ec/{filename}_MAP-TO_mbg-k{{kmer}}-w{{window}}.ms{{minscore}}.clip-ec.stats.tsv.gz',
            filename=ont_read_files
            )
    output:
        summary = 'output/seq_stats/summary/HG00733_ONT_clip-ec.k{kmer}-w{window}.ms{minscore}.stats.tsv',
        dump = 'output/seq_stats/summary/HG00733_ONT_clip-ec.k{kmer}-w{window}.ms{minscore}.stats.pck'
    log:
        'log/output/seq_stats/summary/HG00733_ONT_clip-ec.k{kmer}-w{window}.ms{minscore}.stats.log'
    benchmark:
        'rsrc/output/seq_stats/summary/HG00733_ONT_clip-ec.k{kmer}-w{window}.ms{minscore}.stats.rsrc'
    conda:
        '../../environment/conda/conda_pyscript.yml'
    threads: config['num_cpu_low']
    params:
        script_exec = lambda wildcards: find_script_path('collect_read_stats.py')
    shell:
        '{params.script_exec} --debug --input-files {input.fasta} '
        '--output {output.dump} --summary-output {output.summary} '
        '--num-cpu {threads} --genome-size 3100000000 &> {log}'


rule master:
    input:
        output_files