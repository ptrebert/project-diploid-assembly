localrules: master

DATA_FOLDER = '/beeond/data/hifi'

REFERENCE_FOLDER = '/beeond/data/references'
REFERENCE_ASSEMBLY = 'T2Tv1_T2TC_chm13'

ALIGNMENT_TARGETS = [
    'T2Tv1_T2TC_chm13',
    'T2Tv1_38p13Y_chm13'
]

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


rule count_reference_kmers:
    input:
        fasta = os.path.join(REFERENCE_FOLDER, '{reference}.fasta')
    output:
        kmer_db = directory('output/kmer_db_sample/{reference}.k{kmer_size}.no-hpc.db'),
        rep_kmer = 'output/kmer_db_sample/{reference}.k{kmer_size}.no-hpc.rep-grt09998.txt'
    benchmark:
        'rsrc/output/kmer_db/{reference}.k{kmer_size}.no-hpc.count-dump.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        reference = REFERENCE_ASSEMBLY
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        mem_total_gb = lambda wildcards, attempt: 32 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    params:
        kmer_size = KMER_SIZE,
        zip_threads = config['num_cpu_high']
    shell:
        'meryl count k={params.kmer_size} threads={threads} memory={resources.mem_total_gb} output {output.kmer_db} {input.fasta} && '
        'meryl print greater-than distinct=0.9998 {output.kmer_db} | pigz -p {params.zip_threads} --best > {output.rep_kmer}'


rule count_sequence_kmers:
    input:
        fastq = select_hifi_input
    output:
        kmer_db = directory('output/kmer_db_sample/{sample}.k{kmer_size}.{hpc}.db'),
        rep_kmer = 'output/kmer_db_sample/{sample}.k{kmer_size}.{hpc}.rep-grt09998.txt'
    benchmark:
        'rsrc/output/kmer_db/{sample}.k{kmer_size}.{hpc}.count-dump.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = '(' + '|'.join(MALE_SAMPLES + FEMALE_SAMPLES) + ')'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 262144 * attempt,
        mem_total_gb = lambda wildcards, attempt: 256 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    params:
        kmer_size = KMER_SIZE,
        zip_threads = config['num_cpu_high'],
        use_hpc = lambda wildcards: '' if wildcards.hpc == 'no-hpc' else 'compress'
    shell:
        'meryl count k={params.kmer_size} threads={threads} memory={resources.mem_total_gb} {params.use_hpc} output {output.kmer_db} {input.fastq} && '
        'meryl print greater-than distinct=0.9998 {output.kmer_db} | pigz -p {params.zip_threads} --best > {output.rep_kmer}'


rule build_shared_male_db:
    input:
        kmer_dbs = expand(
            'output/kmer_db_sample/{sample}.k{{kmer_size}}.{{hpc}}.db',
            sample=MALE_SAMPLES
        )
    output:
        kmer_db = directory('output/kmer_db/male-shared.k{kmer_size}.{hpc}.db')
    benchmark:
        'rsrc/output/kmer_db/male-shared.k{kmer_size}.{hpc}.intersect.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_gb = lambda wildcards, attempt: 8 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'meryl intersect-sum threads={threads} memory={resources.mem_total_gb} {input.kmer_dbs} output {output.kmer_db}'


rule build_merged_female_db:
    input:
        kmer_dbs = expand(
            'output/kmer_db_sample/{sample}.k{{kmer_size}}.{{hpc}}.db',
            sample=FEMALE_SAMPLES
        )
    output:
        kmer_db = directory('output/kmer_db/female-merged.k{kmer_size}.{hpc}.db')
    benchmark:
        'rsrc/output/kmer_db/female-merged.k{kmer_size}.{hpc}.union.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 * attempt,
        mem_total_gb = lambda wildcards, attempt: 16 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'meryl union-sum threads={threads} memory={resources.mem_total_gb} {input.kmer_dbs} output {output.kmer_db}'


rule build_male_specific_db:
    input:
        male_kmers = 'output/kmer_db/male-shared.k{kmer_size}.{hpc}.db',
        female_kmers = 'output/kmer_db/female-merged.k{kmer_size}.{hpc}.db',
        reference_kmers = 'output/kmer_db_sample/{}.k{{kmer_size}}.no-hpc.db'.format(REFERENCE_ASSEMBLY),
    output:
        male_specific = directory('output/kmer_db/male-specific.k{kmer_size}.{hpc}.db')
    benchmark:
        'rsrc/output/kmer_db/male-specific.k{kmer_size}.{hpc}.difference.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
        mem_total_gb = lambda wildcards, attempt: 12 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'meryl difference threads={threads} memory={resources.mem_total_gb} '
        '{input.male_kmers} {input.female_kmers} {input.reference_kmers} '
        'output {output.male_specific}'


rule dump_male_specific_kmers:
    input:
        kmer_db = 'output/kmer_db/male-specific.k{kmer_size}.{hpc}.db'
    output:
        txt = 'output/kmer_db/male-specific.k{kmer_size}.{hpc}.txt.gz'
    benchmark:
        'rsrc/output/kmer_db/male-specific.k{kmer_size}.{hpc}.dump.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_total_gb = lambda wildcards, attempt: attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'meryl print {input.kmer_db} | pigz -p {threads} --best > {output.txt}'


rule dump_male_unique_kmers:
    input:
        kmer_db = 'output/kmer_db/male-specific.k{kmer_size}.{hpc}.db'
    output:
        txt = 'output/kmer_db/male-unique.k{kmer_size}.{hpc}.txt.gz'
    benchmark:
        'rsrc/output/kmer_db/male-unique.k{kmer_size}.{hpc}.dump.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_total_gb = lambda wildcards, attempt: attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    params:
        shared_uniq = len(MALE_SAMPLES)
    shell:
        'meryl print [equal-to {params.shared_uniq} {input.kmer_db}] | pigz -p {threads} --best > {output.txt}'


rule extract_ktagged_reads:
    input:
        fastq = select_hifi_input,
        kmer_db = 'output/kmer_db/male-specific.k{kmer_size}.{hpc}.db'
    output:
        ktagged_reads = 'output/ktagged_reads/{sample}.k{kmer_size}.{hpc}.ktagged-reads.fastq.gz',
    benchmark:
        'rsrc/output/ktagged_reads/{sample}.k{kmer_size}.{hpc}.ktagged-reads.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = '(' + '|'.join(MALE_SAMPLES) + ')'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_gb = lambda wildcards, attempt: 8 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    params:
        kmer_size = KMER_SIZE,
        zip_threads = config['num_cpu_high'] - 2,
    shell:
        'meryl-lookup -include -sequence {input.fastq} -output /dev/stdout -mers {input.kmer_db} | '
        'pigz -p {params.zip_threads} --best > {output.ktagged_reads}'


rule count_parental_kmers:
    input:
        fastq = select_hifi_input
    output:
        dump = 'ouput/kmer_dumps/{sample}.k{kmer_size}.yak'
    benchmark:
        'rsrc/ouput/kmer_dumps/{sample}.k{kmer_size}.yak.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_max']
    resources:
        mem_total_mb = lambda wildcards, attempt: 262144 * attempt,
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    params:
        kmer_size = KMER_SIZE,
    shell:
        'yak count -k {params.kmer_size} -b 37 -K 4096m -o {output.dump} {input.fastq}' 


def select_maternal_kmer_dump(wildcards):
    path = 'ouput/kmer_dumps/{sample}.k{kmer_size}.yak'
    mothers = {
        'NA24385': 'NA24143'
    }
    mother = mothers[wildcards.sample]
    formatter = {
        'sample': mother,
        'kmer_size': KMER_SIZE
    }
    return path.format(**formatter)


def select_paternal_kmer_dump(wildcards):
    path = 'ouput/kmer_dumps/{sample}.k{kmer_size}.yak'
    fathers = {
        'NA24385': 'NA24149'
    }
    father = fathers[wildcards.sample]
    formatter = {
        'sample': father,
        'kmer_size': KMER_SIZE
    }
    return path.format(**formatter)


rule compute_hifiasm_trio_assembly:
    input:
        fastq = select_hifi_input,
        mat_yak = select_maternal_kmer_dump,
        pat_yak = select_paternal_kmer_dump
    output:
        assm_dir = directory('output/assemblies/trio_binned/{sample}')
    log:
        'log/output/assemblies/trio_binned/{sample}.hifiasm.log',
    benchmark:
        'rsrc/output/assemblies/trio_binned/{sample}.hifiasm.rsrc',
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_max']
    resources:
        mem_total_mb = lambda wildcards, attempt: 524288 * attempt,
        runtime_hrs = lambda wildcards, attempt: 72 * attempt
    shell:
        'hifiasm -o {output.assm_dir}/{wildcards.sample} -t {threads} -1 {input.pat_yak} -2 {input.mat_yak} {input.fastq} &> {log.hifiasm}'


rule compute_hifiasm_trio_assembly:
    input:
        fastq = select_hifi_input,
    output:
        assm_dir = directory('output/assemblies/non_trio/{sample}')
    log:
        'log/output/assemblies/non_trio/{sample}.hifiasm.log',
    benchmark:
        'rsrc/output/assemblies/non_trio/{sample}.hifiasm.rsrc',
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_max']
    resources:
        mem_total_mb = lambda wildcards, attempt: 524288 * attempt,
        runtime_hrs = lambda wildcards, attempt: 72 * attempt
    shell:
        'hifiasm -o {output.assm_dir}/{wildcards.sample} -t {threads} {input.fastq} &> {log.hifiasm}'


rule master:
    input:
        expand(
            rules.dump_male_specific_kmers.output.txt,
            kmer_size=KMER_SIZE,
            hpc=['is-hpc']
        ),
        expand(
            rules.dump_male_unique_kmers.output.txt,
            kmer_size=KMER_SIZE,
            hpc=['is-hpc']
        ),
        expand(
            rules.extract_ktagged_reads.output.ktagged_reads,
            kmer_size=KMER_SIZE,
            hpc=['is-hpc'],
            sample=MALE_SAMPLES
        ),
        'output/assemblies/trio_binned/NA24385',
        expand(
            'output/assemblies/non_trio/{sample}',
            sample=MALE_SAMPLES
        )

