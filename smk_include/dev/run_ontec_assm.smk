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

output_files = [
    'input/ont/HG00733_EEE_ONT.fa.gz',
    'input/ont/HG00733_1_Guppy_4.2.2_prom.fa.gz',
    'input/ont/HG00733_2_Guppy_4.2.2_prom.fa.gz',
    'input/ont/HG00733_3_Guppy_4.2.2_prom.fa.gz',
    'output/mbg_hifi/HG00733_HiFi.mbg-k5001-w2000.gfa',
    'output/mbg_hifi/HG00733_HiFi.mbg-k2001-w1000.gfa',
    'output/mbg_hifi/HG00733_HiFi.mbg-k501-w100.gfa',
]

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
        mem_per_cpu_mb = lambda wildcards, attempt: int(2048 * attemp // config['num_cpu_low']),
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
        mem_per_cpu_mb = lambda wildcards, attempt: int(2048 * attemp // config['num_cpu_low']),
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt + attempt ** attempt
    params:
        pigz_cpu = lambda wildcards: int(config['num_cpu_low']) - 1
    shell:
        'seqtk seq -S -A -l 120 -C {input.fastq} | pigz -p {params.pigz_cpu} > {output.fasta}'


rule build_hifi_read_graph:
    """
    https://github.com/maickrau/MBG/issues/1

    ### I've used vg: vg view -Fv graph.gfa | vg mod -n -U 100 - | vg view - > blunt-graph.gfa
    """
    input:
        '/gpfs/project/projects/medbioinf/projects/hifi_v13/run_folder/input/fastq/HG00733_hgsvc_pbsq2-ccs_1000.fastq.gz'
    output:
        'output/mbg_hifi/HG00733_HiFi.mbg-k{kmer}-w{window}.gfa'
    log:
        'log/output/mbg_hifi/HG00733_HiFi.mbg-k{kmer}-w{window}.log'
    benchmark:
        'rsrc/output/mbg_hifi/HG00733_HiFi.mbg-k{kmer}-w{window}.rsrc'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(32768 * attemp // config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    shell:
        'MBG -i {input} -o {output} -t {threads} --blunt -k {wildcards.kmer} -w {wildcards.window} &> {log}'


rule master:
    input:
        output_files