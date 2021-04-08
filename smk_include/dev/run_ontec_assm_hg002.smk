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
    'output/seq_stats/input/ont/HG002_giab_ULfastqs_guppy324.stats.tsv.gz',
    'output/seq_stats/input/ont/HG002_ONT_PAD64459_Guppy32.stats.tsv.gz',

]

pattern_reads = 'output/ont_ec/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ev{evalue}.clip-ec.fa.gz'
pattern_graph = 'output/mbg_hifi_clean/HG002_HiFi.mbg-k{kmer}-w{window}.clean.noseq.gfa'
for orf in ont_read_files:
    for k, w in [(5001,2000), (2001,1000), (501,100)]:
        tmp = pattern_reads.format(**{
            'filename': orf,
            'kmer': k,
            'window': w,
            'evalue': '1e-2'
        })
        output_files.append(tmp)

        stats_output = tmp.rsplit('.', 2)[0] + '.stats.tsv.gz'
        stats_output = os.path.join('output', 'seq_stats', stats_output)
        output_files.append(stats_output)

        tmp = pattern_graph.format(**{
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
        filename = '(' + '|'.join(ont_read_files) + ')'
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
        mem_per_cpu_mb = lambda wildcards, attempt: int(122880 + 122880 * attempt // config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 122880 + 122880 * attempt,
        runtime_hrs = lambda wildcards, attempt: 3 * attempt
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
        mem_per_cpu_mb = lambda wildcards, attempt: 16384 * attempt,
        mem_total_mb = lambda wildcards, attempt: 16384 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt ** attempt
    shell:
        'vg view -Fv {input} | vg mod -n -U 100 - 2> {log} | vg view - > {output}'


rule strip_sequences_from_graph:
    input:
        'output/mbg_hifi_clean/HG002_HiFi.mbg-k{kmer}-w{window}.clean.gfa'
    output:
        'output/mbg_hifi_clean/HG002_HiFi.mbg-k{kmer}-w{window}.clean.noseq.gfa'
    benchmark:
        'rsrc/output/mbg_hifi_clean/HG002_HiFi.mbg-k{kmer}-w{window}.noseq.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt ** attempt
    shell:
        'gfatools -S {input} > {output}'



rule ont_error_correction:
    """
    CPU and runtime resources adapted for ultra-slow cluster I/O
    """
    input:
        graph = 'output/mbg_hifi_clean/HG002_HiFi.mbg-k{kmer}-w{window}.clean.gfa',
        reads = 'input/ont/{filename}.fa.gz',
    output:
        gaf = 'output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ev{evalue}.gaf',
        ec_reads = 'output/ont_ec/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ev{evalue}.clip-ec.fa.gz',
    log:
        'log/output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ev{evalue}.ga.log'
    benchmark:
        'rsrc/output/ont_aln/{filename}_MAP-TO_mbg-k{kmer}-w{window}.ev{evalue}.ga.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(94208 * attempt // config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 94208 * attempt,
        runtime_hrs = lambda wildcards, attempt: 72 * attempt
    params:
        evalue = lambda wildcards: float(wildcards.evalue)
    shell:
        'GraphAligner -t {threads} -g {input.graph} -f {input.reads} '
            '-x dbg --E-cutoff {params.evalue} '
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


rule master:
    input:
        output_files