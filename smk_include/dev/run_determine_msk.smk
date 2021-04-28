localrules: master

DATA_FOLDER = '/beeond/data/hifi'

MALE_SAMPLES = [
    'HG00731',
    'HG00512',
    'NA19239',
    'NA24385',
    'NA24149'
]

FEMALE_SAMPLES = [
    'HG00732',
    'HG00733',
    'HG00513',
    'HG00514',
    'NA19238',
    'NA19240',
    'NA24143'
]

KMER_SIZE = 31


def select_hifi_input(wildcards):
    if wildcards.sample in ['NA24143', 'NA24149', 'NA24385']:
        return os.path.join(DATA_FOLDER, '{}_hpg_pbsq2-ccs_1000.fastq.gz'.format(wildcards.sample))
    else:
        return os.path.join(DATA_FOLDER, '{}_hgsvc_pbsq2-ccs_1000.fastq.gz'.format(wildcards.sample))


rule count_sequence_kmers:
    input:
        fastq = select_hifi_input
    output:
        kmer_db = 'output/kmer_db/{sample}.k{kmer_size}.db/'),
        rep_kmer = 'output/kmer_db/{sample}.k{kmer_size}.rep-grt09998.txt')
    benchmark:
        'rsrc/output/kmer_db/{sample}.k{kmer_size}.count-dump.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 262144 * attempt,
        mem_total_gb = lambda wildcards, attempt: 256 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    params:
        kmer_size = KMER_SIZE,
        zip_threads = config['num_cpu_high']
    shell:
        'meryl count k={params.kmer_size} threads={threads} memory={resources.mem_total_gb} output {output.kmer_db} {input.fastq} && '
        'meryl print greater-than distinct=0.9998 {output.kmer_db} | pigz -p {params.zip_threads} --best > {output.rep_kmer}'


rule build_shared_male_db:
    input:
        kmer_dbs = expand(
            'output/kmer_db/{sample}.k{{kmer_size}}.db/'),
            sample=MALE_SAMPLES
        )
    output:
        kmer_db = 'output/kmer_db/male-shared.k{kmer_size}.db/'
    benchmark:
        'rsrc/output/kmer_db/male-shared.k{kmer_size}.intersect.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 262144 * attempt,
        mem_total_gb = lambda wildcards, attempt: 256 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    shell:
        'meryl intersect-sum threads={threads} memory={resources.mem_total_gb} {input.kmer_dbs} output {output.kmer_db}'


rule build_merged_female_db:
    input:
        kmer_dbs = expand(
            'output/kmer_db/{sample}.k{{kmer_size}}.db/'),
            sample=FEMALE_SAMPLES
        )
    output:
        kmer_db = 'output/kmer_db/female-merged.k{kmer_size}.db/'
    benchmark:
        'rsrc/output/kmer_db/female-merged.k{kmer_size}.union.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 262144 * attempt,
        mem_total_gb = lambda wildcards, attempt: 256 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    shell:
        'meryl union-sum threads={threads} memory={resources.mem_total_gb} {input.kmer_dbs} output {output.kmer_db}'


rule build_male_specific_db:
    input:
        male_kmers = 'output/kmer_db/male-shared.k{kmer_size}.db/',
        female_kmers = 'output/kmer_db/female-merged.k{kmer_size}.db/'
    output:
        male_specific = 'output/kmer_db/male-specific.k{kmer_size}.db/'
    benchmark:
        'rsrc/output/kmer_db/male-specific.k{kmer_size}.difference.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 262144 * attempt,
        mem_total_gb = lambda wildcards, attempt: 256 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    shell:
        'meryl difference threads={threads} memory={resources.mem_total_gb} {input.male_kmers} {input.female_kmers} output {output.male_specific}'


rule dump_male_specific_kmers:
    input:
        kmer_db = 'output/kmer_db/male-specific.k{kmer_size}.db/'
    output:
        txt = 'output/kmer_db/male-specific.k{kmer_size}.txt.gz'
    benchmark:
        'rsrc/output/kmer_db/male-specific.k{kmer_size}.dump.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        mem_total_gb = lambda wildcards, attempt: 32 * attempt,
        runtime_hrs = lambda wildcards, attempt: 24 * attempt
    shell:
        'meryl print threads={threads} memory={resources.mem_total_gb} {input.kmer_db} | pigz -p {threads} --best > {output.txt}'


rule master:
    input:
        expand(
            rules.dump_male_specific_kmers.output.txt,
            kmer_size=KMER_SIZE
        )
