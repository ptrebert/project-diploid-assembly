localrules: master

### INPUTFILES

ont_read_files = [
    'HG002_giab_ULfastqs_guppy324',
    'HG002_ONT_PAD64459_Guppy32',
]

output_files = [
    'input/ont/HG002_giab_ULfastqs_guppy324.fa.gz',
    'input/ont/HG002_ONT_PAD64459_Guppy32.fa.gz',
    'output/mbg_hifi/HG002_HiFi.mbg-k5001-w2000.gfa',
    'output/mbg_hifi/HG002_HiFi.mbg-k2001-w1000.gfa',
    'output/mbg_hifi/HG002_HiFi.mbg-k501-w100.gfa',
]

pattern = 'output/ont_ec/{filename}_MAP-TO_mbg-k{kmer}-w{window}.clip-ec.fa.gz'
for orf in ont_read_files:
    for k, w in [(5001,2000), (2001,1000), (501,100)]:
        tmp = pattern.format(**{
            'filename': orf,
            'kmer': k,
            'window': w
        })
        output_files.append(tmp)


rule clean_hpg_ont:
    input:
        fastq = '/beeond/data/ont/{filename}.fastq.gz'
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
        '/beeond/projects/lansdorp/run_folder/input/fastq/NA24385_hpg_pbsq2-ccs_1000.fastq.gz'
    output:
        'output/mbg_hifi/HG002_HiFi.mbg-k{kmer}-w{window}.gfa'
    log:
        'log/output/mbg_hifi/HG002_HiFi.mbg-k{kmer}-w{window}.log'
    benchmark:
        'rsrc/output/mbg_hifi/HG002_HiFi.mbg-k{kmer}-w{window}.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(91136 * attempt // config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 91136 * attempt,
        runtime_hrs = lambda wildcards, attempt: 12 * attempt
    shell:
        'MBG -i {input} -o {output} -t {threads} --blunt -k {wildcards.kmer} -w {wildcards.window} &> {log}'


rule clean_mbg_graph:
    """
    https://github.com/maickrau/MBG/issues/1

    ### I've used vg: vg view -Fv graph.gfa | vg mod -n -U 100 - | vg view - > blunt-graph.gfa
    """
    input:
        'output/mbg_hifi/HG002_HiFi.mbg-k{kmer}-w{window}.gfa'
    output:
        'output/mbg_hifi_clean/HG002_HiFi.mbg-k{kmer}-w{window}.clean.gfa'
    log:
        'log/output/mbg_hifi_clean/HG002_HiFi.mbg-k{kmer}-w{window}.clean.log'
    benchmark:
        'rsrc/output/mbg_hifi_clean/HG002_HiFi.mbg-k{kmer}-w{window}.clean.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt ** attempt
    shell:
        'vg view -Fv {input} | vg mod -n -U 100 - 2> {log} | vg view - > {output}'


rule ont_error_correction:
    input:
        graph = 'output/mbg_hifi_clean/HG002_HiFi.mbg-k{kmer}-w{window}.clean.gfa',
        reads = 'input/ont/{filename}.fa.gz',
    output:
        gaf = 'output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.gaf',
        ec_reads = 'output/ont_ec/{filename}_MAP-TO_mbg-k{kmer}-w{window}.clip-ec.fa.gz',
    log:
        'log/output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ga.log'
    benchmark:
        'rsrc/output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ga.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(91136 * attempt // config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 91136 * attempt,
        runtime_hrs = lambda wildcards, attempt: 16 * attempt
    shell:
        'GraphAligner -t {threads} -g {input.graph} -f {input.reads} -x dbg --corrected-clipped-out {output.ec_reads} -a {output.gaf} &> {log}'


rule master:
    input:
        output_files