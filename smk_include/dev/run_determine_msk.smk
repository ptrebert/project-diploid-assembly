localrules: create_ceny_url_files, master

REFERENCE_FOLDER = '/gpfs/project/projects/medbioinf/data/references'
REFERENCE_ASSEMBLY = 'T2Tv11_T2TC_chm13'

ASSEMBLY_INPUT = '/gpfs/project/projects/medbioinf/data/share/globus/sig_chrY/assemblies/freeze_v1/graphs'

ALIGNMENT_TARGETS = [
    'T2Tv11_T2TC_chm13',
    'T2Tv11_38p13Y_chm13',
    'T2Tv11_38p13Ycen_chm13'
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
    'NA24143',
    'HG03125',
    'HG02818',
    'HG03486',
    'NA12878'
]

HIFI_SAMPLES = MALE_SAMPLES + FEMALE_SAMPLES

WMAP_KMER_LONG_READS = 15
# deprecated: use minimap2 (v2.20) for contig-to-reference mapping
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

wildcard_constraints:
    hpc = '(ishpc|nohpc)'


def select_assembly_input(wildcards):

    infix = _select_hifi_input(wildcards)
    infix = os.path.basename(infix).rsplit('.', 2)[0]

    tig_map = {
        'rutg': 'r_utg',
        'pctg': 'p_ctg',
        'actg': 'a_ctg'
    }
    tigs = tig_map.get(wildcards.tigs, None)
    if tigs is None:
        raise ValueError('Cannot process tig value: {}'.format(str(wildcards)))

    assembly_graph = os.path.join(
        ASSEMBLY_INPUT,
        infix,
        infix + '.{}.gfa'.format(tigs)
    )
    return assembly_graph


def _select_hifi_input(wildcards):

    input_folder = 'input/hifi_reads_filtered'
    if wildcards.sample in ['NA24143', 'NA24149', 'NA24385']:
        readset = '{}_hpg_pbsq2-ccs_1000.fastq.gz'.format(wildcards.sample)
    elif wildcards.sample in ['NA12878']:
        readset = '{}_giab_pbsq2-ccs_1000.fastq.gz'.format(wildcards.sample)
    else:
        readset = '{}_hgsvc_pbsq2-ccs_1000.fastq.gz'.format(wildcards.sample)
    return os.path.join(input_folder, readset)


def _select_reference_input(wildcards):

    if not wildcards.sample in ALIGNMENT_TARGETS:
        raise ValueError('Presumably invalid reference: {}'.format(wildcards.sample))
    return os.path.join(REFERENCE_FOLDER, wildcards.sample + '.fasta')


def _select_read_input(wildcards):

    if wildcards.sample == 'HG02982A0':
        seq_files =  ['references/assemblies/HG02982A0.fasta']
    elif wildcards.sample == 'RP11CENYILL':
        seq_files = expand(
            'references/downloads/reads/RP11CENY_{platform}_{accession}_{mate}.fastq.gz',
            platform=['ILL'],
            accession=RP11CENY_ILL_ACCESSIONS,
            mate=[1, 2]
        )
    elif wildcards.sample == 'RP11CENYONT':
        seq_files = expand(
            'references/downloads/reads/RP11CENY_{platform}_{accession}_{mate}.fastq.gz',
            platform=['ONT'],
            accession=RP11CENY_ONT_ACCESSIONS,
            mate=[1]
        )
    else:
        raise ValueError(str(wildcards))
    return seq_files


def select_sequence_input(wildcards):

    if 'T2T' in wildcards.sample:
        return _select_reference_input(wildcards)
    elif wildcards.sample in ['HG02982A0', 'RP11CENYILL', 'RP11CENYONT']:
        return _select_read_input(wilkdcards)
    else:
        return _select_hifi_input(wildcards)


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


rule run_hic_assemblies:
    input:
        'output/assemblies/layout/hifi_hic/HG00731.done',
        'output/assemblies/layout/hifi_hic/HG002.done'


rule run_all:
    input:
        # align male reference reads/sequences to sample assembly graphs
        # (if alignable, prefer raw unitg graph)
        expand(
            'output/alignments/reads_to_graph/{male_reads}_MAP-TO_{sample}.{tigs}.gaf',
            male_reads=['GRCh38_chrY', 'HG02982_A0'],  # 'RP11CENY_ONT' --- atm, seems impossible to ga these reads
            sample=MALE_SAMPLES,
            tigs=['rutg', 'pctg', 'actg']
        ),
        # dump sample assembly graph tigs to BED and align tigs as reads to reference
        expand(
            'output/alignments/tigs_to_reference/{sample}.{tigs}_MAP-TO_{reference}.{content}.bed',
            sample=MALE_SAMPLES,
            tigs=['rutg'],
            reference=['T2Tv11_38p13Ycen_chm13'],
            content=['filt', 'unmap']
        ),
        # identify male-specific long reads and align to reference
        expand(
            'output/alignments/ktagged_to_ref/{sample}.k{msk_kmer}.{hpc}.ktg-{share_type}_MAP-TO_{reference}.wmap-k{wmap_kmer}.{content}.bed',
            sample=MALE_SAMPLES,
            msk_kmer=[KMER_SIZE],
            wmap_kmer=[WMAP_KMER_LONG_READS],
            reference=['T2Tv11_38p13Ycen_chm13'],
            hpc=['ishpc'],
            content=['filt', 'unmap'],
            share_type=['specific']  # 'specific-uniq', 'union-uniq': "unique" k-mer sets hardly contain any info (indeed empty for union-uniq)
        ),


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
        'references/downloads/reads/RP11CENY_{platform}_{accession}_{mate}.fastq.gz'
    shell:
        'wget --no-verbose -O {output} -i {input}'
    

rule merge_rp11ceny_ont_reads:
    input:
        expand(
            'references/downloads/reads/RP11CENY_ONT_{accession}_1.fastq.gz',
            accession=RP11CENY_ONT_ACCESSIONS
        )
    output:
        'references/reads/RP11CENY_ONT.reads.fasta'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    shell:
        'pigz -d -c -p {threads} {input} | seqtk seq -C -A > {output}'


rule load_ceny_assembly:
    output:
        'references/downloads/Jain2018_MF741337.fasta'
    shell:
        'wget --no-verbose -O {output} "https://www.ebi.ac.uk/ena/browser/api/fasta/MF741337.1?download=true" '


rule preprocess_rp11ceny_assembly:
    input:
        fasta = 'references/downloads/Jain2018_MF741337.fasta'
    output:
        fasta = 'references/assemblies/MF741337.fasta'
    run:
        import io
        out_buffer = io.StringIO()
        with open(input.fasta, 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    out_buffer.write('>chrYCEN\n')
                else:
                    out_buffer.write(line.strip())
        
        with open(output.fasta, 'w') as dump:
            _ = dump.write(out_buffer.getvalue())


rule extend_male_t2t_assembly:
    input:
        cen = 'references/assemblies/MF741337.fasta',
        ref = ancient(os.path.join(REFERENCE_FOLDER, 'T2Tv11_38p13Y_chm13.fasta'))
    output:
        fasta = os.path.join(REFERENCE_FOLDER, 'T2Tv11_38p13Ycen_chm13.fasta')
    shell:
        'cat {input.ref} {input.cen} > {output.fasta}'


def count_kmer_runtime(wildcards, attempt):

    if wildcards.sample in HIFI_SAMPLES:
        return 24 * attempt
    elif 'ONT' in wildcards.sample:
        return attempt * attempt * attempt
    else:
        return attempt * attempt


def count_kmer_memory(wildcards, attempt, unit='mb'):

    if wildcards.sample in HIFI_SAMPLES:
        mem = 262144
    elif 'ONT' in wildcards.sample:
        mem = 73728
    elif 'T2T' in wildcards.sample:
        mem = 32768
    else:
        mem = 4096
    if unit == 'gb':
        mem = int(mem / 1024)
    return mem * attempt


rule count_kmers:
    input:
        sequence = select_sequence_input
    output:
        kmer_db = directory('output/kmer_db_sample/{sample}.k{kmer_size}.{hpc}.db'),
        rep_kmer = 'output/kmer_db_sample/{sample}.k{kmer_size}.{hpc}.rep-grt09998.txt'
    benchmark:
        'rsrc/output/kmer_db/{sample}.k{kmer_size}.{hpc}.count-dump.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        kmer_size = '(' + '|'.join([str(WMAP_KMER_LONG_READS), str(WMAP_KMER_ASSM_CTG), str(KMER_SIZE)]) + ')'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: count_kmer_memory(wildcards, attempt),
        mem_total_gb = lambda wildcards, attempt: count_kmer_memory(wildcards, attempt, 'gb'),
        runtime_hrs = lambda wildcards, attempt: count_kmer_runtime(wildcards, attempt)
    params:
        kmer_size = KMER_SIZE,
        use_hpc = lambda wildcards: '' if wildcards.hpc == 'nohpc' else 'compress'
    shell:
        'meryl count k={threads} threads={threads} memory={resources.mem_total_gb} {params.use_hpc} output {output.kmer_db} {input.sequence} && '
        'meryl print greater-than distinct=0.9998 {output.kmer_db} > {output.rep_kmer}'


rule build_rp11ceny_specific_db:
    input:
        ont_db = 'output/kmer_db_sample/RP11CENYONT.k{kmer_size}.is-hpc.db',
        ill_db = 'output/kmer_db_sample/RP11CENYILL.k{kmer_size}.no-hpc.db',
        female_db = 'output/kmer_db/female-merged.k{kmer_size}.is-hpc.db',
        ref_db = 'output/kmer_db_sample/{}.k{{kmer_size}}.is-hpc.db'.format(REFERENCE_ASSEMBLY),
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
        a0_db = 'output/kmer_db_sample/HG02982A0.k{kmer_size}.{hpc}.db',
        female_db = 'output/kmer_db/female-merged.k{kmer_size}.{hpc}.db',
        ref_db = 'output/kmer_db_sample/{}.k{{kmer_size}}.{{hpc}}.db'.format(REFERENCE_ASSEMBLY),
    output:
        kmer_db = directory('output/kmer_db/HG02982A0-specific.k{kmer_size}.{hpc}.db')
    benchmark:
        'rsrc/output/kmer_db/HG02982A0-specific.k{kmer_size}.{hpc}.build.rsrc'
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
        reference_kmers = 'output/kmer_db_sample/{}.k{{kmer_size}}.{{hpc}}.db'.format(REFERENCE_ASSEMBLY),
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
        a0_kmers = 'output/kmer_db/HG02982A0-specific.k{kmer_size}.{hpc}.db',
        rp11_kmers = 'output/kmer_db/RP11CENY-specific.k{kmer_size}.{hpc}.db'
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
        a0_kmers = 'output/kmer_db/HG02982A0-specific.k{kmer_size}.{hpc}.db',
        rp11_kmers = 'output/kmer_db/RP11CENY-specific.k{kmer_size}.{hpc}.db'
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
        fastq = select_sequence_input,
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


rule ktagged_reads_to_linear_reference_alignment:
    input:
        reads = 'output/ktagged_reads/{sample}.k{msk_kmer}.{hpc}.ktg-{share_type}.fastq.gz',
        reference = ancient(os.path.join(REFERENCE_FOLDER, '{reference}.fasta')),
        ref_repkmer = 'output/kmer_db_sample/{reference}.k{wmap_kmer}.{hpc}.rep-grt09998.txt',
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
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 24576 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt ** 5,
        mem_sort_mb = 4096,
        align_threads = config['num_cpu_high'] - config['num_cpu_low'],
        sort_threads = config['num_cpu_low'],
    params:
        individual = lambda wildcards: wildcards.sample,
        readgroup_id = lambda wildcards: '{}_k{}_{}_{}'.format(
            wildcards.sample,
            wildcards.msk_kmer,
            wildcards.hpc,
            wildcards.share_type.replace('-', '')
        ),
        tempdir = lambda wildcards: os.path.join(
            'temp', 'winnowmap', 'ktagref',
            wildcards.reference,
            wildcards.sample,
            wildcards.hpc,
            wildcards.share_type,
            'MSK' + str(wildcards.msk_kmer),
            'WMAP' + str(wildcards.wmap_kmer)
        )
    shell:
        'rm -rfd {params.tempdir} ; mkdir -p {params.tempdir} && '
        'winnowmap -W {input.ref_repkmer} -k {wildcards.wmap_kmer} -t {resources.align_threads} -a -x map-pb  '
        '-R "@RG\\tID:{params.readgroup_id}\\tSM:{params.individual}" --secondary=no '
        '{input.reference} {input.reads} | '
        'samtools sort -m {resources.mem_sort_mb}M -@ {resources.sort_threads} -T {params.tempdir} -O BAM > {output.bam} ; '
        'rm -rfd {params.tempdir}'


rule graph_align_male_reference_reads:
    input:
        reads = 'references/reads/{male_reads}.reads.fasta',
        graph = select_assembly_input
    output:
        gaf = 'output/alignments/reads_to_graph/{male_reads}_MAP-TO_{sample}.{tigs}.gaf'
    log:
        'log/output/alignments/reads_to_graph/{male_reads}_MAP-TO_{sample}.{tigs}.ga.log'
    benchmark:
        'rsrc/output/alignments/reads_to_graph/{male_reads}_MAP-TO_{sample}.{tigs}.ga.rsrc'
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


rule convert_tigs_gfa_to_fasta:
    input:
        gfa = select_assembly_input
    output:
        fasta = 'output/assemblies/{sample}.{tigs}.fasta',
        rc_map = 'output/assemblies/{sample}.{tigs}.read-contig.map',
        stats = 'output/assemblies/{sample}.{tigs}.contig.stats',
    log:
        'log/output/assemblies/{sample}.{tigs}.gfa-convert.log'
    benchmark:
        'run/output/assemblies/{sample}.{tigs}.gfa-convert' + '.t{}.rsrc'.format(config['num_cpu_low'])
    conda:
        '../../environment/conda/conda_pyscript.yml'
    wildcard_constraints:
        sample = '(' + '|'.join(MALE_SAMPLES) + ')',
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script_exec = lambda wildcards: find_script_path('gfa_to_fasta.py')
    shell:
        '{params.script_exec} --gfa {input.gfa} --n-cpus {threads} '
        '--out-fasta {output.fasta} --out-map {output.rc_map} --out-stats {output.stats} &> {log}'


rule contig_to_linear_reference_alignment:
    """
    Important that this rule relies on the unimap backports
    for better contig-to-ref alignment available in minimap2 2.20+
    """
    input:
        contigs = 'output/assemblies/{sample}.{tigs}.fasta',
        ref_fasta = ancient(os.path.join(REFERENCE_FOLDER, '{reference}.fasta'))
    output:
        'output/alignments/tigs_to_reference/{sample}.{tigs}_MAP-TO_{reference}.psort.raw.bam'
    log:
        mm = 'log/output/alignments/tigs_to_reference/{sample}.{tigs}_MAP-TO_{reference}.psort.raw.mm.log',
        st = 'log/output/alignments/tigs_to_reference/{sample}.{tigs}_MAP-TO_{reference}.psort.raw.st.log',
    benchmark:
        'rsrc/output/alignments/tigs_to_reference/{sample}.{tigs}_MAP-TO_{reference}.mmap.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = '(' + '|'.join(MALE_SAMPLES) + ')'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((32768 + 32768 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: 24 * attempt,
        align_threads = config['num_cpu_high'] - config['num_cpu_low'],
        mem_sort_mb = 4096,
        sort_threads = config['num_cpu_low']
    params:
        individual = lambda wildcards: wildcards.sample,
        readgroup_id = lambda wildcards: '{}_{}'.format(wildcards.sample, wildcards.tigs),
        tempdir = lambda wildcards: os.path.join(
            'temp', 'minimap', 'ctgref',
            wildcards.reference,
            wildcards.sample,
            wildcards.tigs
        ),
    shell:
        'rm -rfd {params.tempdir} ; mkdir -p {params.tempdir} && '
        'minimap2 -t {resources.align_threads} '
            '--secondary=no --eqx -Y -ax asm5 '
            '-R "@RG\\tID:{params.readgroup_id}\\tSM:{params.individual}" '
            '{input.ref_fasta} {input.contigs} 2> {log.mm} | '
            'samtools sort -@ {resources.sort_threads} -m {resources.mem_sort_mb}M -T {params.tempdir} -O BAM > {output} 2> {log.st} ; '
        'rm -rfd {params.tempdir}'


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
    df = df.loc[df['chrom'].isin(['chrX', 'chrY', 'chrYCEN']), ['chrom', 'read_name', 'length']].copy()
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
            'output/alignments/ktagged_to_ref/{{sample}}.k{{msk_kmer}}.{{hpc}}.{ktag_set}_MAP-TO_{reference}.wmap-k15.{aln_type}.bed',
            ktag_set=['ktg-specific', 'ktg-specific-uniq', 'ktg-union'],
            reference=['T2Tv11_38p13Ycen_chm13'],
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

        info_columns = ('chrY', 'chrYCEN', 'chrX', 'chrUn', 'specific', 'specific_uniq', 'union')
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
        'output/ktagged_reads/gonosomal/{sample}.k{msk_kmer}.{hpc}.{assembly}.subset-table.tsv',
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
            subset = tig_map.loc[tig_map['read'].isin(read_table.loc[selector, 'read']), :].copy()
            subset['label'] = label
            subset['color'] = color
            subset['overlap'] = 1
            subset = subset[['contig', 'label', 'color', 'overlap']]
            subsets.append(subset)
        
        subsets = pd.concat(subsets, axis=0)
        subsets.drop_duplicates(keep='first', inplace=True)
        with open(output[0], 'w') as table:
            subsets.to_csv(table, sep='\t', header=False, index=False)


rule strip_sequence_from_graph:
    input:
        '{file_path}/{file_name}.gfa'
    output:
        '{file_path}/{file_name}.noseq.gfa'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'gfatools view -S {input} > {output}'


def tag_graph_alignments(gaf_path):

    translation_table = dict((i,i) for i in '1234567890')
    translation_table['>'] = ' '
    translation_table['<'] = ' '
    translation_table = str.maketrans(translation_table)

    tagged_tigs = set()
    path_members = set()

    with open(gaf_path, 'r') as gaf:
        for line in gaf:
            tig = line.split()[5]
            tig = set(tig.translate(translation_table).split())
            tagged_tigs = tagged_tigs.union(tig)
            if len(tig) > 1:
                path_members = path_members.union(tig)
    return tagged_tigs, path_members

    
def tag_linear_alignments(bed_path):

    tagged_tigs = set()
    cen_tigs = set()
    with open(bed_path, 'r') as bed:
        for line in bed:
            if not line.startswith('chrY'):
                continue
            parts = line.split()
            chrom = parts[0]
            tig = parts[3]
            if chrom == 'chrYCEN':
                cen_tigs.add(tig)
            else:
                tagged_tigs.add(tig)
    return tagged_tigs, cen_tigs
    

rule tag_raw_unitigs_by_alignment:
    input:
        ref_reads = 'output/alignments/reads_to_graph/GRCh38_chrY_MAP-TO_{sample}.{tigs}.gaf',
        a0_reads = 'output/alignments/reads_to_graph/HG02982_A0_MAP-TO_{sample}.{tigs}.gaf',
        lin_aln = 'output/alignments/tigs_to_reference/{sample}.{tigs}_MAP-TO_T2Tv11_38p13Ycen_chm13.filt.bed',
        #cen_reads = 'output/alignments/reads_to_graph/{male_reads}_MAP-TO_{sample}.rutg.gaf'
    output:
        'output/tagged_{tigs}/{sample}_{tigs}.ga_hg38_A0.ln_hg38_cen.tsv'
    run:
        import pandas as pd
        import numpy as np

        ref_tags, ref_paths = tag_graph_alignments(input.ref_reads)
        a0_tags, a0_paths = tag_graph_alignments(input.a0_reads)
        lin_tags, cen_tags = tag_linear_alignments(input.lin_aln)

        labels = [
            'ga_hg38_tag', 'ga_hg38_path',
            'ga_A0_tag', 'ga_A0_path',
            'ln_hg38', 'ln_RP11CEN'
        ]

        all_tagged = [ref_tags, ref_paths, a0_tags, a0_paths, lin_tags, cen_tags]

        merged_tags = set().union(*all_tagged)
        num_rows = len(merged_tags)
        num_columns = len(labels)
        df = pd.DataFrame(
            np.zeros((num_rows, num_columns), dtype=np.int8),
            index=pd.Index(sorted(merged_tags), name=wildcards.tigs),
            columns=labels
        )

        for tigs, label in zip(all_tagged, labels):
            df.loc[df.index.isin(tigs), label] = 1

        df.to_csv(output[0], sep='\t', index=True, header=True)


rule create_simple_table_gfasubset:
    input:
        'output/tagged_{tigs}/{sample}_{tigs}.ga_hg38_A0.ln_hg38_cen.tsv'
    output:
        'output/tagged_{tigs}/{sample}_{tigs}.ga_hg38_A0.ln_hg38_cen.simple.tsv'
    run:
        import pandas as pd

        df = pd.read_csv(input[0], sep='\t', header=0, index_col=None)
        df['reg_label'] = 'unset'
        df['reg_color'] = 'unset'
        df['overlap_bp'] = 1

        select_cen = (df['ln_RP11CEN'] == 1)
        df.loc[select_cen, 'reg_label'] = 'CEN'
        df.loc[select_cen, 'reg_color'] = 'blue'

        select_empty = df['reg_label'] == 'unset'
        select_a0 = (df['ga_A0_tag'] == 1) | (df['ga_A0_path'] == 1)
        select_a0 = select_a0 & select_empty
        df.loc[select_a0, 'reg_label'] = 'A0'
        df.loc[select_a0, 'reg_color'] = 'red'

        select_empty = df['reg_label'] == 'unset'
        df.loc[select_empty, 'reg_label'] = 'hg38'
        df.loc[select_empty, 'reg_color'] = 'black'

        assert (df['reg_label'] != 'unset').all(), 'Region missed'
        df[[wildcards.tigs, 'reg_label', 'reg_color', 'overlap_bp']].to_csv(
            output[0],
            sep='\t',
            header=False,
            index=False
        )


rule subset_assembly_graph:
    input:
        table = 'output/tagged_{tigs}/{sample}_{tigs}.ga_hg38_A0.ln_hg38_cen.simple.tsv',
        graph = select_assembly_input
    output:
        graph = 'output/subset_graphs/{sample}.{tigs}.chrY.gfa',
        table = 'output/subset_graphs/{sample}.{tigs}.chrY.csv',
    log:
        'log/output/subset_graphs/{sample}.{tigs}.chrY.subset.log',
    conda:
        '../../environment/conda/conda_pyscript.yml'
    params:
        script_exec = lambda wildcards: find_script_path('gfa_subset.py'),
        comp_cov = 10
    shell:
        '{params.script_exec} --debug --input-gfa {input.graph} --simple-table --input-table {input.table} '
        '--output-gfa {output.graph} --output-table {output.table} --component-tag-coverage {params.comp_cov} &> {log}'


rule hifiasm_hic_assembly_hg002:
    """
    """
    input:
        hifi_reads = '/gpfs/project/projects/medbioinf/data/share/globus/sig_chrY/clean_hifi_reads/NA24385_hpg_pbsq2-ccs_1000.fastq.gz',
        hic_reads1 = '/gpfs/project/projects/medbioinf/data/share/globus/sig_chrY/HiC/HG002/HG002_comb_R1_1.fastq.gz',
        hic_reads2 = '/gpfs/project/projects/medbioinf/data/share/globus/sig_chrY/HiC/HG002/HG002_comb_R2_2.fastq.gz',
    output:
        assm_done = 'output/assemblies/layout/hifi_hic/HG002.done'
    log:
        hifiasm = 'log/output/assemblies/layout/hifi_hic/HG002.hifiasm.log'
    benchmark:
        'rsrc/output/assemblies/layout/hifi_hic/HG002.hifiasm.t24.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'       
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 + 83968 * attempt,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    params:
        prefix = lambda wildcards, output: output.assm_done.rsplit('.', 1)[0],
    shell:
        'hifiasm -o {params.prefix} -t {threads} --h1 {input.hic_reads1} --h2 {input.hic_reads2} {input.hifi_reads} &> {log.hifiasm} '
        ' && touch {output.assm_done}'


rule hifiasm_hic_assembly_hg00731:
    """
    """
    input:
        hifi_reads = '/gpfs/project/projects/medbioinf/data/share/globus/sig_chrY/clean_hifi_reads/HG00731_hgsvc_pbsq2-ccs_1000.fastq.gz',
        hic_reads1 = '/gpfs/project/projects/medbioinf/data/share/globus/sig_chrY/HiC/HG00731/HG00731_biosamples_1-2_comb_1.fastq.gz',
        hic_reads2 = '/gpfs/project/projects/medbioinf/data/share/globus/sig_chrY/HiC/HG00731/HG00731_biosamples_1-2_comb_2.fastq.gz',
    output:
        assm_done = 'output/assemblies/layout/hifi_hic/HG00731.done'
    log:
        hifiasm = 'log/output/assemblies/layout/hifi_hic/HG00731.hifiasm.log'
    benchmark:
        'rsrc/output/assemblies/layout/hifi_hic/HG00731.hifiasm.t24.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'       
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 + 83968 * attempt,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    params:
        prefix = lambda wildcards, output: output.assm_done.rsplit('.', 1)[0],
    shell:
        'hifiasm -o {params.prefix} -t {threads} --h1 {input.hic_reads1} --h2 {input.hic_reads2} {input.hifi_reads} &> {log.hifiasm} '
        ' && touch {output.assm_done}'
