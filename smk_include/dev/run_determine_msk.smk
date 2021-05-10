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


rule determine_ktagged_overlapping_reads:
    input:
        fastq = select_hifi_input,
        ktagged = 'output/ktagged_reads/{sample}.k{kmer_size}.{hpc}.ktagged-reads.fastq.gz',
    output:
        read_ovl = 'output/alignments/ktag_to_hifi/{sample}.k{kmer_size}.{hpc}.hifi-ovl.paf',
    log:
        'log/output/alignments/ktag_to_hifi/{sample}.k{kmer_size}.{hpc}.hifi-ovl.log'
    benchmark:
        'rsrc/output/alignments/ktag_to_hifi/{sample}.k{kmer_size}.{hpc}.hifi-ovl.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = '(' + '|'.join(MALE_SAMPLES) + ')'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: 36 * attempt
    shell:
        'minimap2 -H -x ava-pb -X -o {output.read_ovl} -t {threads} {input.fastq} {input.ktagged}'


rule count_parental_kmers:
    input:
        fastq = select_hifi_input
    output:
        dump = 'output/kmer_dumps/{sample}.k{kmer_size}.yak'
    benchmark:
        'rsrc/output/kmer_dumps/{sample}.k{kmer_size}.yak.rsrc'
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
    path = 'output/kmer_dumps/{sample}.k{kmer_size}.yak'
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
    path = 'output/kmer_dumps/{sample}.k{kmer_size}.yak'
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
    """
    v0.15.1-r328
    """
    input:
        fastq = select_hifi_input,
        mat_yak = select_maternal_kmer_dump,
        pat_yak = select_paternal_kmer_dump
    output:
        assm_out = multiext(
            'output/assemblies/trio_binned/{sample}/{sample}.dip',
            '.hap1.p_ctg.gfa', '.hap1.p_ctg.lowQ.bed', '.hap1.p_ctg.noseq.gfa',
            '.hap2.p_ctg.gfa', '.hap2.p_ctg.lowQ.bed', '.hap2.p_ctg.noseq.gfa',
            '.p_utg.gfa', '.p_utg.lowQ.bed', '.p_utg.noseq.gfa',
            '.r_utg.gfa', '.r_utg.lowQ.bed', '.r_utg.noseq.gfa'
        )
    log:
        'log/output/assemblies/trio_binned/{sample}.hifiasm.log',
    benchmark:
        'rsrc/output/assemblies/trio_binned/{sample}.hifiasm.rsrc',
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_max']
    resources:
        mem_total_mb = lambda wildcards, attempt: 114688 + 114688 * attempt,
        runtime_hrs = lambda wildcards, attempt: 24 + 8 * attempt
    shell:
        'hifiasm -o {output.assm_dir}/{wildcards.sample} -t {threads} -1 {input.pat_yak} -2 {input.mat_yak} {input.fastq} &> {log}'


rule compute_hifiasm_nontrio_assembly:
    """
    v0.15.1-r328
    """
    input:
        fastq = select_hifi_input,
    output:
        assm_out = multiext(
            'output/assemblies/non_trio/{sample}/{sample}.bp',
            '.hap1.p_ctg.gfa', '.hap1.p_ctg.lowQ.bed', '.hap1.p_ctg.noseq.gfa',
            '.hap2.p_ctg.gfa', '.hap2.p_ctg.lowQ.bed', '.hap2.p_ctg.noseq.gfa',
            '.p_ctg.gfa', '.p_ctg.lowQ.bed', '.p_ctg.noseq.gfa',
            '.p_utg.gfa', '.p_utg.lowQ.bed', '.p_utg.noseq.gfa',
            '.r_utg.gfa', '.r_utg.lowQ.bed', '.r_utg.noseq.gfa'
        )
    log:
        'log/output/assemblies/non_trio/{sample}.hifiasm.log',
    benchmark:
        'rsrc/output/assemblies/non_trio/{sample}.hifiasm.rsrc',
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_max']
    resources:
        mem_total_mb = lambda wildcards, attempt: 73728 + 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 24 + 8 * attempt
    params:
        prefix = lambda wildcards: os.path.join(
            'output/assemblies/non_trio', wildcards.sample, wildcards.sample)
    shell:
        'hifiasm -o {params.prefix} -t {threads} {input.fastq} &> {log}'


rule align_male_reference_reads:
    input:
        reads = 'references/reads/{male_reads}.reads.fasta',
        graph = 'output/assemblies/{assm_mode}/{sample}/{sample}.{mode_id}.{hap}.p_ctg.gfa'
    output:
        gaf = 'output/alignments/reads_to_graph/{male_reads}_MAP-TO_{sample}_{assm_mode}.{mode_id}.{hap}.p_ctg.gaf'
    log:
        'log/output/alignments/reads_to_graph/{male_reads}_MAP-TO_{sample}_{assm_mode}.{mode_id}.{hap}.p_ctg.ga.log'
    benchmark:
        'rsrc/output/alignments/reads_to_graph/{male_reads}_MAP-TO_{sample}_{assm_mode}.{mode_id}.{hap}.p_ctg.ga.rsrc'
    wildcard_constraints:
        sample = '(' + '|'.join(MALE_SAMPLES) + ')'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_max']  # change to high on HILBERT
    resources:
        mem_total_mb = lambda wildcards, attempt: 65536 + 65536 * attempt,
        runtime_hrs = lambda wildcards, attempt: 6 + 6 * attempt
    shell:
        'GraphAligner --verbose -x vg -t {threads} '
        '--multimap-score-fraction 0.99 --min-alignment-score 1000 '
        '-g {input.graph} -f {input.reads} -a {output.gaf} &> {log}'


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
        ),
        expand(
            'output/alignments/reads_to_graph/{male_reads}_MAP-TO_{sample}_{assm_mode}.dip.{hap}.p_ctg.gaf',
            assm_mode=['trio_binned'],
            male_reads=['GRCh38_chrY', 'HG02982_A0'],
            sample=['NA24385'],
            hap=['hap1', 'hap2']
        ),
        expand(
            'output/alignments/reads_to_graph/{male_reads}_MAP-TO_{sample}_{assm_mode}.bp.{hap}.p_ctg.gaf',
            assm=['non_trio']
            male_reads=['GRCh38_chrY', 'HG02982_A0'],
            sample=MALE_SAMPLES,
            hap=['hap1', 'hap2']
        ),
        expand(
            'output/alignments/ktag_to_hifi/{sample}.k{kmer_size}.{hpc}.hifi-ovl.paf',
            kmer_size=[KMER_SIZE],
            sample=MALE_SAMPLES,
            hpc=['is-hpc']
        )
