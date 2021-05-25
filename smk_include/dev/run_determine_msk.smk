localrules: mock_index_reads, create_ceny_url_files, master

DATA_FOLDER = '/beeond/data/hifi'

REFERENCE_FOLDER = '/beeond/data/references'
REFERENCE_ASSEMBLY = 'T2Tv11_T2TC_chm13'

ALIGNMENT_TARGETS = [
    'T2Tv1_T2TC_chm13',
    'T2Tv1_38p13Y_chm13',
    'T2Tv11_T2TC_chm13',
    'T2Tv11_38p13Y_chm13'
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

WMAP_KMER_LONG_READS = 15
# deprecated: use unimap for contig-to-reference mapping
WMAP_KMER_ASSM_CTG = 19

# k-mer size for male-specific k-mer search
KMER_SIZE = 31


RP11CENY_ONT_ACCESSIONS = [
    'SRR5902337',
    'SRR5902338',
    'SRR5902339',
    'SRR5902340',
    'SRR5902341',
    'SRR5902342',
    'SRR5902343',
    'SRR5902344',
    'SRR5902345',
    'SRR5902346',
]

RP11CENY_ILL_ACCESSIONS = [
    'SRR5902347',
    'SRR5902348',
    'SRR5902349',
    'SRR5902350',
    'SRR5902351',
    'SRR5902352',
    'SRR5902353',
    'SRR5902354',
    'SRR5902355'
]



def select_hifi_input(wildcards):
    if wildcards.sample in ['NA24143', 'NA24149', 'NA24385']:
        return os.path.join(DATA_FOLDER, '{}_hpg_pbsq2-ccs_1000.fastq.gz'.format(wildcards.sample))
    else:
        return os.path.join(DATA_FOLDER, '{}_hgsvc_pbsq2-ccs_1000.fastq.gz'.format(wildcards.sample))


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


rule load_haplogroupA_assembly:
    """
    ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/ulg/ULGL01.fasta.gz
    """
    output:
        'references/downloads/HG02982A0.ULGL01.PRJEB28143.fasta.gz'
    shell:
        'wget --no-verbose -O {output} ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/ulg/ULGL01.fasta.gz'


rule preprocess_haplogroupA_assembly:
    input:
        'references/downloads/HG02982A0.ULGL01.PRJEB28143.fasta.gz'
    output:
        fasta = 'references/assemblies/HG02982A0.fasta',
        stats = 'references/assemblies/HG02982A0.contig.stats',
    run:
        import gzip
        import collections as col
        import operator as op

        seq_buffer = col.defaultdict(str)

        current_contig = None
        with gzip.open(input[0], 'rt') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    current_contig = line.strip().strip('>')
                else:
                    seq_buffer[current_contig] += line.strip().upper()
        
        seq_stats = []
        for contig, ctg_seq in seq_buffer.items():
            seq_stats.append((len(ctg_seq), contig, col.Counter(ctg_seq)))
        
        with open(output.fasta, 'w') as fasta:
            for _, contig, stats in sorted(seq_stats, reverse=True):
                simple_contig = contig.split()[-1]
                simple_contig = simple_contig.split('_')[0]
                assert simple_contig.startswith('tig00'), 'Renaming failed: {}'.format(simple_contig)
                fasta_header = '>HG02982_A0_' + simple_contig
                _ = fasta.write('{}\n{}\n'.format(fasta_header, seq_buffer[contig]))
        
        get_nuc_counts = op.itemgetter(*('A', 'C', 'G', 'T'))
        with open(output.stats, 'w') as table:
            _ = table.write('\t'.join(['name', 'length', 'A', 'C', 'G', 'T']) + '\n')
            for length, contig, stats in sorted(seq_stats, reverse=True):
                simple_contig = contig.split()[-1]
                simple_contig = simple_contig.split('_')[0]
                assert simple_contig.startswith('tig00'), 'Renaming failed: {}'.format(simple_contig)
                a,c,g,t = get_nuc_counts(stats)
                _ = table.write('\t'.join(list(map(str, [simple_contig, length, a, c, g, t]))) + '\n')
        # END OF RUN BLOCK


rule create_ceny_url_files:
    input:
        'annotation/filereport_read_run_PRJNA397218_json.txt'
    output:
        'annotation/RP11CENY_{platform}_{accession}_{mate}.url'
    run:
        import json

        with open(input[0], 'r') as dump:
            metadata = json.load(dump)

        for md in metadata:
            if md['run_accession'] != wildcards.accession:
                continue
            remote_urls = md['fastq_ftp'].split(';')
            check_mate_num = '_' + str(wildcards.mate) + '.fastq.gz'
            for r in remote_urls:
                if not r.endswith(check_mate_num):
                    continue
                full_url = 'ftp://' + r
                with open(output[0], 'w') as dump:
                    _ = dump.write(full_url)


rule load_rp11ceny_reads:
    input:
        'annotation/RP11CENY_{platform}_{accession}_{mate}.url'
    output:
        'output/reads/RP11CENY_{platform}_{accession}_{mate}.fastq.gz'
    shell:
        'wget --no-verbose -O {output} -i {input}'
    

rule merge_rp11ceny_ont_reads:
    input:
        expand(
            'output/reads/RP11CENY_ONT_{accession}_1.fastq.gz',
            accession=RP11CENY_ONT_ACCESSIONS
        )
    output:
        'references/reads/RP11CENY_ONT.reads.fasta'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    shell:
        'pigz -d -c -p {threads} {input} | seqtk seq -C -A > {output}'


rule mock_index_reads:
    input:
        os.path.join(REFERENCE_FOLDER, '{reference}.fasta')
    output:
        temp(os.path.join(REFERENCE_FOLDER, 'mock_index', '{reference}.index_read.fasta'))
    run:
        with open(output[0], 'w') as fasta:
            _ = fasta.write('>index_read\n')
            _ = fasta.write('ACGTACGT\n')


rule create_unimap_index:
    """
    NB: index compatibility (k-mer size default: 21)
    """
    input:
        ref = os.path.join(REFERENCE_FOLDER, '{reference}.fasta'),
        reads = os.path.join(REFERENCE_FOLDER, 'mock_index', '{reference}.index_read.fasta')
    output:
        umi = os.path.join(REFERENCE_FOLDER, '{reference}.umi'),
    benchmark:
        'run/references/indexing/{reference}.umi.rsrc',
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: 2
    resources:
        runtime_hrs = lambda wildcards, attempt: max(0, attempt - 1),
        mem_total_mb = lambda wildcards, attempt: 16384 + 16384 * attempt
    shell:
        'unimap -d {output} -x asm20 -t {threads} -o /dev/null {input.ref} {input.reads}'


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
        reference = '(' + '|'.join(ALIGNMENT_TARGETS) + ')',
        kmer_size = '(' + '|'.join([str(WMAP_KMER_LONG_READS), str(WMAP_KMER_ASSM_CTG), str(KMER_SIZE)]) + ')'
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


def select_read_file_input(wildcards):

    if wildcards.sample == 'HG02982A0':
        seq_files =  ['references/assemblies/HG02982A0.fasta']
    elif wildcards.sample == 'RP11CENYILL':
        seq_files = expand(
            'output/reads/RP11CENY_{platform}_{accession}_{mate}.fastq.gz',
            platform=['ILL'],
            accession=RP11CENY_ILL_ACCESSIONS,
            mate=[1, 2]
        )
    elif wildcards.sample == 'RP11CENYONT':
        seq_files = expand(
            'output/reads/RP11CENY_{platform}_{accession}_{mate}.fastq.gz',
            platform=['ONT'],
            accession=RP11CENY_ONT_ACCESSIONS,
            mate=[1]
        )
    else:
        raise ValueError(str(wildcards))
    return seq_files


rule count_read_kmers:
    input:
        seq_files = select_read_file_input
    output:
        kmer_db = directory('output/kmer_db_sample/{sample}.k{kmer_size}.{hpc}.db'),
        rep_kmer = 'output/kmer_db_sample/{sample}.k{kmer_size}.{hpc}.rep-grt09998.txt'
    benchmark:
        'rsrc/output/kmer_db/{sample}.k{kmer_size}.{hpc}.count-dump.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = '(RP11CENYONT|RP11CENYILL|HG02982A0)'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 if attempt < 2 else 73728,
        mem_total_gb = lambda wildcards, attempt: 2 if attempt < 2 else 72,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        kmer_size = KMER_SIZE,
        zip_threads = config['num_cpu_high'],
        use_hpc = lambda wildcards: '' if wildcards.hpc == 'no-hpc' else 'compress'
    shell:
        'meryl count k={params.kmer_size} threads={threads} memory={resources.mem_total_gb} {params.use_hpc} output {output.kmer_db} {input.seq_files} && '
        'meryl print greater-than distinct=0.9998 {output.kmer_db} | pigz -p {params.zip_threads} --best > {output.rep_kmer}'


rule build_rp11ceny_specific_db:
    input:
        ont_db = 'output/kmer_db_sample/RP11CENYONT.k{kmer_size}.is-hpc.db',
        ill_db = 'output/kmer_db_sample/RP11CENYILL.k{kmer_size}.no-hpc.db',
        female_db = 'output/kmer_db/female-merged.k{kmer_size}.is-hpc.db',
        ref_db = 'output/kmer_db_sample/{}.k{{kmer_size}}.no-hpc.db'.format(REFERENCE_ASSEMBLY),
    output:
        kmer_db = directory('output/kmer_db/RP11CENY-specific.k{kmer_size}.is-hpc.db')
    benchmark:
        'rsrc/output/kmer_db/RP11CENY-specific.k{kmer_size}.is-hpc.build.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_gb = lambda wildcards, attempt: 8 * attempt,
        runtime_hrs = lambda wildcards, attempt: 8 * attempt
    shell:
        'meryl threads={threads} memory={resources.mem_total_gb} difference [intersect-sum {input.ont_db} {input.ill_db}] '
            ' {input.female_db} {input.ref_db} output {output.kmer_db}'


rule build_a0_specific_db:
    input:
        a0_db = 'output/kmer_db_sample/HG02982A0.k{kmer_size}.is-hpc.db',
        female_db = 'output/kmer_db/female-merged.k{kmer_size}.is-hpc.db',
        ref_db = 'output/kmer_db_sample/{}.k{{kmer_size}}.no-hpc.db'.format(REFERENCE_ASSEMBLY),
    output:
        kmer_db = directory('output/kmer_db/HG02982A0-specific.k{kmer_size}.is-hpc.db')
    benchmark:
        'rsrc/output/kmer_db/HG02982A0-specific.k{kmer_size}.is-hpc.build.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_gb = lambda wildcards, attempt: 8 * attempt,
        runtime_hrs = lambda wildcards, attempt: 8 * attempt
    shell:
        'meryl threads={threads} memory={resources.mem_total_gb} difference {input.a0_db} '
            ' {input.female_db} {input.ref_db} output {output.kmer_db}'


rule build_shared_male_db:
    input:
        kmer_dbs = expand(
            'output/kmer_db_sample/{sample}.k{{kmer_size}}.{{hpc}}.db',
            sample=MALE_SAMPLES
        )
    output:
        kmer_db = directory('output/kmer_db/male-shared.k{kmer_size}.{hpc}.db')
    benchmark:
        'rsrc/output/kmer_db/male-shared.k{kmer_size}.{hpc}.build.rsrc'
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
        'rsrc/output/kmer_db/female-merged.k{kmer_size}.{hpc}.build.rsrc'
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
        'rsrc/output/kmer_db/male-specific.k{kmer_size}.{hpc}.build.rsrc'
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


rule build_male_specific_uniq_db:
    """
    k-mers occurring once in all male HiFi data sets, but never
    in female HiFi reads or in the T2T reference
    """
    input:
        male_kmers = 'output/kmer_db/male-specific.k{kmer_size}.{hpc}.db',
    output:
        male_specific_uniq = directory('output/kmer_db/male-specific-uniq.k{kmer_size}.{hpc}.db')
    benchmark:
        'rsrc/output/kmer_db/male-specific-uniq.k{kmer_size}.{hpc}.build.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
        mem_total_gb = lambda wildcards, attempt: 12 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        male_uniq = len(MALE_SAMPLES),
    shell:
        'meryl threads={threads} memory={resources.mem_total_gb} '
        'equal-to {params.male_uniq} {input.male_kmers} '
        'output {output.male_specific_uniq}'


rule build_male_union_db:
    input:
        male_kmers = 'output/kmer_db/male-specific.k{kmer_size}.{hpc}.db',
        a0_kmers = 'output/kmer_db/HG02982A0-specific.k{kmer_size}.is-hpc.db',
        rp11_kmers = 'output/kmer_db/RP11CENY-specific.k{kmer_size}.is-hpc.db'
    output:
        male_union = directory('output/kmer_db/male-union.k{kmer_size}.{hpc}.db')
    benchmark:
        'rsrc/output/kmer_db/male-union.k{kmer_size}.{hpc}.build.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
        mem_total_gb = lambda wildcards, attempt: 12 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'meryl union-sum threads={threads} memory={resources.mem_total_gb} '
        '{input.male_kmers} {input.a0_kmers} {input.rp11_kmers} '
        'output {output.male_union}'


rule build_male_union_uniq_db:
    """
    Because of the bias introduced by RP11-CENY k-mers, this should
    predominantly contain k-mers corresponding to centromere sequence
    """
    input:
        male_kmers = 'output/kmer_db/male-specific.k{kmer_size}.{hpc}.db',
        a0_kmers = 'output/kmer_db/HG02982A0-specific.k{kmer_size}.is-hpc.db',
        rp11_kmers = 'output/kmer_db/RP11CENY-specific.k{kmer_size}.is-hpc.db'
    output:
        male_union_uniq = directory('output/kmer_db/male-union-uniq.k{kmer_size}.{hpc}.db')
    benchmark:
        'rsrc/output/kmer_db/male-union.k{kmer_size}.{hpc}.build.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
        mem_total_gb = lambda wildcards, attempt: 12 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        male_uniq = len(MALE_SAMPLES),
        union_uniq = len(MALE_SAMPLES) + 2
    shell:
        'meryl threads={threads} memory={resources.mem_total_gb} '
        'intersect-sum '
        '[equal-to {params.male_uniq} {input.male_kmers}] '
        '[equal-to 1 {input.a0_kmers}] '
        '[equal-to 1 {input.rp11_kmers}] '
        'output {output.male_union_uniq}'


rule dump_male_kmers:
    input:
        kmer_db = 'output/kmer_db/male-{share_type}.k{kmer_size}.{hpc}.db'
    output:
        txt = 'output/kmer_db/male-{share_type}.k{kmer_size}.{hpc}.txt.gz'
    benchmark:
        'rsrc/output/kmer_db/male-{share_type}.k{kmer_size}.{hpc}.dump-txt.rsrc'
    wildcard_constraints:
        share_type = '(specific|union|specific-uniq|union-uniq)'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_total_gb = lambda wildcards, attempt: attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'meryl print {input.kmer_db} | pigz -p {threads} --best > {output.txt}'


rule extract_ktagged_reads:
    input:
        fastq = select_hifi_input,
        kmer_db = 'output/kmer_db/male-{share_type}.k{kmer_size}.{hpc}.db'
    output:
        ktagged_reads = 'output/ktagged_reads/{sample}.k{kmer_size}.{hpc}.ktg-{share_type}.fastq.gz',
    benchmark:
        'rsrc/output/ktagged_reads/{sample}.k{kmer_size}.{hpc}.ktg-{share_type}.rsrc'
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


rule align_ktagged_reads_to_reference:
    input:
        reads = 'output/ktagged_reads/{sample}.k{msk_kmer}.{hpc}.ktg-{share_type}.fastq.gz',
        reference = os.path.join(REFERENCE_FOLDER, '{reference}.fasta'),
        ref_repkmer = 'output/kmer_db_sample/{reference}.k{wmap_kmer}.no-hpc.rep-grt09998.txt',
    output:
        bam = 'output/alignments/ktagged_to_ref/{sample}.k{msk_kmer}.{hpc}.ktg-{share_type}_MAP-TO_{reference}.wmap-k{wmap_kmer}.psort.raw.bam'
    log:
        'log/output/alignments/ktagged_to_ref/{sample}.k{msk_kmer}.{hpc}.ktg-{share_type}_MAP-TO_{reference}.wmap-k{wmap_kmer}.log'
    benchmark:
        'rsrc/output/alignments/ktagged_to_ref/{sample}.k{msk_kmer}.{hpc}.ktg-{share_type}_MAP-TO_{reference}.wmap-k{wmap_kmer}.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = '(' + '|'.join(MALE_SAMPLES) + ')'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt ** 5,
        mem_sort_mb = 4096,
    params:
        individual = lambda wildcards: wildcards.sample,
        readgroup_id = lambda wildcards: '{}_k{}_{}_{}'.format(
            wildcards.sample,
            wildcards.msk_kmer,
            wildcards.hpc.replace('-', ''),
            wildcards.share_type.replace('-', '')
        ),
        align_threads = config['num_cpu_medium'] - 4,
        sort_threads = 4
    shell:
        'winnowmap -W {input.ref_repkmer} -k {wildcards.wmap_kmer} -t {params.align_threads} -a -x map-pb  '
        '-R "@RG\\tID:{params.readgroup_id}\\tSM:{params.individual}" --secondary=no '
        '{input.reference} {input.reads} | '
        'samtools sort -m {resources.mem_sort_mb}M -@ {params.sort_threads} -O BAM > {output.bam}'


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
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 65536 + 65536 * attempt,
        runtime_hrs = lambda wildcards, attempt: 36 * attempt if 'ONT' in wildcards.male_reads else 6 + 6 * attempt
    shell:
        'GraphAligner --verbose -x vg -t {threads} '
        '--multimap-score-fraction 1 --min-alignment-score 500 '
        '-g {input.graph} -f {input.reads} -a {output.gaf} &> {log}'


rule convert_primary_contigs_gfa_to_fasta:
    input:
        gfa = 'output/assemblies/{assm_mode}/{sample}/{sample}.{mode_id}.{hap}.p_ctg.gfa',
    output:
        fasta = 'output/assemblies/{sample}.{assm_mode}.{mode_id}.{hap}.p_ctg.fasta',
        rc_map = 'output/assemblies/{sample}.{assm_mode}.{mode_id}.{hap}.p_ctg.read-contig.map',
        stats = 'output/assemblies/{sample}.{assm_mode}.{mode_id}.{hap}.p_ctg.contig.stats',
    log:
        'log/output/assemblies/{sample}.{assm_mode}.{mode_id}.{hap}.p_ctg.gfa-convert.log'
    benchmark:
        'run/output/assemblies/{sample}.{assm_mode}.{mode_id}.{hap}.p_ctg.gfa-convert' + '.t{}.rsrc'.format(config['num_cpu_low'])
    conda:
        '../../environment/conda/conda_pyscript.yml'
    wildcard_constraints:
        sample = '(' + '|'.join(MALE_SAMPLES) + ')',
        hap = '(hap1|hap2)'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = lambda wildcards: find_script_path('gfa_to_fasta.py')
    shell:
        '{params.script_exec} --gfa {input.gfa} --n-cpus {threads} '
        '--out-fasta {output.fasta} --out-map {output.rc_map} --out-stats {output.stats}'


rule convert_raw_unitigs_gfa_to_fasta:
    input:
        gfa = 'output/assemblies/{assm_mode}/{sample}/{sample}.{mode_id}.r_utg.gfa',
    output:
        fasta = 'output/assemblies/{sample}.{assm_mode}.{mode_id}.r_utg.fasta',
        rc_map = 'output/assemblies/{sample}.{assm_mode}.{mode_id}.r_utg.read-contig.map',
        stats = 'output/assemblies/{sample}.{assm_mode}.{mode_id}.r_utg.contig.stats',
    log:
        'log/output/assemblies/{sample}.{assm_mode}.{mode_id}.r_utg.gfa-convert.log'
    benchmark:
        'run/output/assemblies/{sample}.{assm_mode}.{mode_id}.r_utg.gfa-convert' + '.t{}.rsrc'.format(config['num_cpu_low'])
    conda:
        '../../environment/conda/conda_pyscript.yml'
    wildcard_constraints:
        sample = '(' + '|'.join(MALE_SAMPLES) + ')'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = lambda wildcards: find_script_path('gfa_to_fasta.py')
    shell:
        '{params.script_exec} --gfa {input.gfa} --n-cpus {threads} '
        '--out-fasta {output.fasta} --out-map {output.rc_map} --out-stats {output.stats}'


rule unimap_contig_to_known_reference_alignment:
    input:
        contigs = 'output/assemblies/{sample}.{assembly}.fasta',
        ref_fa = os.path.join(REFERENCE_FOLDER, '{reference}.umi')
    output:
        'output/alignments/tigs_to_reference/{sample}.{assembly}_MAP-TO_{reference}.psort.raw.bam'
    benchmark:
        'rsrc/output/alignments/tigs_to_reference/{sample}.{assembly}_MAP-TO_{reference}.umap.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = '(' + '|'.join(MALE_SAMPLES) + ')'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((16384 + 32768 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 16384 + 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: 36 * attempt,
        mem_sort_mb = 8192
    params:
        individual = lambda wildcards: wildcards.sample,
        readgroup_id = lambda wildcards: '{}_{}'.format(wildcards.sample, wildcards.assembly),
        tempdir = lambda wildcards: os.path.join(
            'temp', 'unimap', wildcards.sample,
            wildcards.assembly
        )
    shell:
        'rm -rfd {params.tempdir} ; mkdir -p {params.tempdir} && '
        'unimap -t {threads} -2 '
            '--secondary=no --eqx -Y -ax asm5 '
            '-R "@RG\\tID:{params.readgroup_id}\\tSM:{params.individual}" '
            '{input.ref_fa} {input.contigs} | '
            'samtools sort -m {resources.mem_sort_mb}M -T {params.tempdir} -O BAM > {output}'


rule dump_sequence_to_reference_alignment_to_bed:
    input:
        'output/alignments/{aln_path}/{aln_file}.psort.raw.bam'
    output:
        'output/alignments/{aln_path}/{aln_file}.filt.bed'
    benchmark:
        'rsrc/output/alignments/{aln_path}/{aln_file}.bed-filter.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt
    params:
        # not (QC fail | dup | unmapped | secondary)
        discard_flag = 1796  # includes supp. alignments
    shell:
        'samtools view -u -q 1 -F {params.discard_flag} -@ {threads} {input} | bedtools bamtobed -i /dev/stdin > {output}'


rule dump_unmapped_sequence_to_bed:
    input:
        'output/alignments/{aln_path}/{aln_file}.psort.raw.bam'
    output:
        'output/alignments/{aln_path}/{aln_file}.unmap.bed'
    benchmark:
        'rsrc/output/alignments/{aln_path}/{aln_file}.bed-unmap.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt
    params:
        # not (QC fail | dup | unmapped | secondary)
        select_flag = 4
    shell:
        'samtools view -u -f {params.select_flag} -@ {threads} {input} | bedtools bamtobed -i /dev/stdin > {output}'


def _read_mapped_ktag_dump(file_path):
    import pandas as pd

    df = pd.read_csv(
        file_path,
        sep='\t',
        header=None,
        names=['chrom', 'start', 'end', 'read_name', 'mapq', 'orientation']
    )
    df['length'] = (df['end'] - df['start']).astype(int)
    df = df.loc[df['chrom'].isin(['chrX', 'chrY']), ['chrom', 'read_name', 'length']].copy()
    return df


def _read_unmapped_ktag_dump(file_path):

    with open(file_path, 'r') as table:
        for line in table:
            if line.strip():
                raise ValueError('non-empty unmapped dump')
    return []


rule select_gonosomal_reads:
    input:
        bed_files = expand(
            'output/alignments/ktagged_to_ref/{{sample}}.k{{msk_kmer}}.{{hpc}}.{ktag_set}_MAP-TO_T2Tv11_38p13Y_chm13.wmap-k15.{aln_type}.bed',
            ktag_set=['ktg-specific', 'ktg-specific-uniq', 'ktg-union'],
            aln_type=['filt', 'unmap']
        )
    output:
        reads_table = 'output/ktagged_reads/gonosomal/{sample}.k{msk_kmer}.{hpc}.tsv',
    run:
        import pandas as pd
        import collections as col
        import operator as op

        read_buffer = col.defaultdict(col.Counter)

        for bed_file in input.bed_files:
            ktag_set = bed_file.split('_MAP-TO_')[0].split('.')[-1]
            ktag_set = ktag_set.split('-', 1)[1].replace('-', '_')
            if bed_file.endswith('.filt.bed'):
                read_info = _read_mapped_ktag_dump(bed_file)
                for chrom, read_name, aln_length in read_info.itertuples(index=False):
                    read_buffer[read_name][ktag_set] += 1
                    read_buffer[read_name][chrom] += aln_length
            elif bed_file.endswith('.unmap.bed'):
                read_info = _read_unmapped_ktag_dump(bed_file)
                for r in read_info:
                    read_buffer[r][ktag_set] += 1
                    read_buffer[r]['chrUn'] += 1
            else:
                raise ValueError

        info_columns = ('chrY', 'chrX', 'chrUn', 'specific', 'specific_uniq', 'union')
        get_read_info = op.itemgetter(*info_columns)
        line_template = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
        with open(output.reads_table, 'w') as table:
            _ = table.write(line_template.format('read', *info_columns))
            for read, read_infos in read_buffer.items():
                _ = table.write(line_template.format(read, *get_read_info(read_infos)))


rule build_gonosome_gfa_subset_table:
    input:
        reads_table = 'output/ktagged_reads/gonosomal/{sample}.k{msk_kmer}.{hpc}.tsv',
        read_tig_map = 'output/assemblies/{sample}.{assembly}.read-contig.map'
    output:
        'output/ktagged_reads/gonosomal/{sample}.k{msk_kmer}.{hpc}.subset-table.tsv',
    run:
        import pandas as pd
        import numpy as np

        read_table = pd.read_csv(input.reads_table, sep='\t', header=0)
        tig_map = pd.read_csv(input.read_tig_map, sep='\t', header=0)

        select_x = np.logical_and(read_table['chrX'] > 0, read_table['chrY'] == 0)
        select_y = np.logical_and(read_table['chrX'] == 0, read_table['chrY'] > 0)
        select_xy = np.logical_and(read_table['chrX'] > 0, read_table['chrY'] > 0)
        # ignore unmapped for now

        subset_labels = ['X', 'Y', 'XY']
        subset_colors = ['pink', 'blue', 'black']
        subset_selectors = [
            select_x,
            select_y,
            select_xy
        ]

        subsets = []

        for selector, label, color in zip(subset_selectors, subset_labels, subset_colors):
            subset = tig_map.loc[tig['read'].isin(read_table.loc[selector, 'read']), :].copy()
            subset['label'] = label
            subset['color'] = color
            subset['overlap'] = 1
            subset = subset[['contig', 'label', 'color']]
            subsets.append(subset)
        
        subsets = pd.concat(subsets, axis=0)
        with open(output[0], 'w') as table:
            subsets.to_csv(table, sep='\t', header=False)


rule master:
    input:
        expand(
            rules.dump_male_kmers.output.txt,
            kmer_size=KMER_SIZE,
            hpc=['is-hpc'],
            share_type=['specific', 'union', 'specific-uniq', 'union-uniq']
        ),
        expand(
            rules.extract_ktagged_reads.output.ktagged_reads,
            kmer_size=KMER_SIZE,
            hpc=['is-hpc'],
            sample=MALE_SAMPLES,
            share_type=['specific', 'union', 'specific-uniq', 'union-uniq']
        ),
        'output/assemblies/trio_binned/NA24385',
        expand(
            'output/assemblies/non_trio/{sample}',
            sample=MALE_SAMPLES
        ),
        expand(
            'output/alignments/reads_to_graph/{male_reads}_MAP-TO_{sample}_{assm_mode}.dip.{hap}.p_ctg.gaf',
            assm_mode=['trio_binned'],
            male_reads=['GRCh38_chrY', 'HG02982_A0', 'RP11CENY_ONT'],
            sample=['NA24385'],
            hap=['hap1', 'hap2']
        ),
        expand(
            'output/alignments/reads_to_graph/{male_reads}_MAP-TO_{sample}_{assm_mode}.bp.{hap}.p_ctg.gaf',
            assm_mode=['non_trio'],
            male_reads=['GRCh38_chrY', 'HG02982_A0', 'RP11CENY_ONT'],
            sample=MALE_SAMPLES,
            hap=['hap1', 'hap2']
        ),
        # expand(
        #     'output/alignments/ktag_to_hifi/{sample}.k{kmer_size}.{hpc}.hifi-ovl.paf',
        #     kmer_size=[KMER_SIZE],
        #     sample=MALE_SAMPLES,
        #     hpc=['is-hpc']
        # ),
        expand(
            'output/alignments/tigs_to_reference/{sample}.{assembly}_MAP-TO_{reference}.{content}.bed',
            sample=MALE_SAMPLES,
            assembly=[
                'non_trio.bp.hap1.p_ctg',
                'non_trio.bp.hap2.p_ctg',
                'non_trio.bp.r_utg',
            ],
            reference=['T2Tv11_38p13Y_chm13'],
            content=['filt', 'unmap']
        ),
        expand(
            'output/alignments/tigs_to_reference/{sample}.{assembly}_MAP-TO_{reference}.{content}.bed',
            sample=['NA24385'],
            assembly=[
                'trio_binned.dip.hap1.p_ctg',
                'trio_binned.dip.hap2.p_ctg',
                'trio_binned.dip.r_utg',
            ],
            reference=['T2Tv11_38p13Y_chm13'],
            content=['filt', 'unmap']
        ),
        expand(
            'output/alignments/ktagged_to_ref/{sample}.k{msk_kmer}.{hpc}.ktg-{share_type}_MAP-TO_{reference}.wmap-k{wmap_kmer}.{content}.bed',
            sample=MALE_SAMPLES,
            msk_kmer=[KMER_SIZE],
            wmap_kmer=[WMAP_KMER_LONG_READS],
            reference=['T2Tv11_38p13Y_chm13'],
            hpc=['is-hpc'],
            content=['filt', 'unmap'],
            share_type=['specific', 'union']  # 'specific-uniq', 'union-uniq': "unique" k-mer sets hardly contain any info (indeed empty for union-uniq)
        ),