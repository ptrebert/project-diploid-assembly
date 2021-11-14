

rule hifiasm_xypar_targeted_assembly:
    input:
        fastq = expand(
            'output/references/{sample_long}.XYPAR.reads.fastq.gz',
            sample_long=[
                'AFR-YRI-Y117-M_NA19239',
                'AMR-PUR-PR05-M_HG00731',
                'EAS-CHS-SH032-M_HG00512',
                'EUR-ASK-3140-M_NA24385'   
            ]
        )
    output:
        primary_unitigs = 'output/target_assembly/xypar_reads/xypar.p_utg.gfa',
        primary_contigs = 'output/target_assembly/xypar_reads/xypar.p_ctg.gfa',
        raw_unitigs = 'output/target_assembly/xypar_reads/xypar.r_utg.gfa',
        discard = multiext(
                'output/target_assembly/xypar_reads/xypar',
                '.a_ctg.gfa', '.a_ctg.noseq.gfa',
                '.p_ctg.noseq.gfa', '.p_utg.noseq.gfa',
                '.r_utg.noseq.gfa'
            ),
        ec_reads = 'output/target_assembly/xypar_reads/xypar.ec.fa',
        reads_ava = 'output/target_assembly/xypar_reads/xypar.ovlp.paf'
    log:
        hifiasm = 'log/output/target_assembly/xypar_reads.hifiasm.log',
    benchmark:
        os.path.join('rsrc/output/target_assembly/xypar_reads.hifiasm' + '.t{}.rsrc'.format(config['num_cpu_medium']))
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        prefix = lambda wildcards, output: output.primary_contigs.rsplit('.', 2)[0],
    shell:
        'hifiasm -o {params.prefix} -t {threads} --write-ec --write-paf --primary {input.fastq} &> {log.hifiasm}'


def select_chry_reads(wildcards):
    """
    Why does this exist?
    - avoid carrying wildcard for input sequence type (FASTA vs FASTQ)
    """
    if wildcards.read_type == 'HIFIAF':
        seq_type = 'fastq'
    else:
        seq_type = 'fasta'
    template = 'output/read_subsets/chry/{sample_info}_{sample}_{read_type}.chrY-reads.{mapq}.{seq_type}.gz'
    formatter = dict(wildcards)
    formatter['seq_type'] = seq_type
    reads_path = template.format(**formatter)
    return reads_path


rule hifiasm_chry_targeted_assembly:
    input:
        reads = select_chry_reads
    output:
        primary_unitigs = 'output/target_assembly/chry_reads/hifiasm/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.p_utg.gfa',
        raw_unitigs = 'output/target_assembly/chry_reads/hifiasm/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.r_utg.gfa',
    log:
        hifiasm = 'log/output/target_assembly/chry_reads/{sample_info}_{sample}_{read_type}.{mapq}.hifiasm.log',
    benchmark:
        'rsrc/output/target_assembly/chry_reads/{sample_info}_{sample}_{read_type}.{mapq}.hifiasm.t24.rsrc',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = CONSTRAINT_ALL_SAMPLES
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 16384 + 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 ** attempt
    params:
        prefix = lambda wildcards, output: output.primary_unitigs.rsplit('.', 2)[0],
        purge_stringency = 3,
        n_hap = 1,
        #genome_size = '60m' --- with hifiasm 0.16.1, this can be set, but since chrYq is quite
        # variable in size, leave it to the default "auto" and hope for the best
    shell:
        'hifiasm -o {params.prefix} -t {threads} --write-paf --primary '
            '-l {params.purge_stringency} --n-hap {params.n_hap} '
            '{input.reads} &> {log.hifiasm}'


rule mbg_chry_targeted_assembly:
    """
    sif = ancient('mbg.master.sif'),
    sif = ancient('mbg.UnitigResolve.sif'),
    'module load Singularity && singularity exec '
        '--bind /:/hilbert {input.sif} '
    """
    input:
        reads = select_chry_reads
    output:
        graph = 'output/target_assembly/chry_reads/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.k{kmer}-w{window}-r{resolvek}.gfa',
        paths = 'output/target_assembly/chry_reads/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.k{kmer}-w{window}-r{resolvek}.gaf',
    log:
        'log/output/target_assembly/chry_reads/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.k{kmer}-w{window}-r{resolvek}.log',
    benchmark:
        'rsrc/output/target_assembly/chry_reads/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.k{kmer}-w{window}-r{resolvek}.t24.rsrc',
    wildcard_constraints:
        read_type = 'HIFIEC',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt,
    shell:
        'MBG -i {input.reads} -t {threads} '
            '-k {wildcards.kmer} -w {wildcards.window} --resolve-maxk {wildcards.resolvek} '
            '--error-masking no '
            '--out {output.graph} --output-sequence-paths {output.paths} &> {log}'


rule lja_chry_targeted_assembly:
    input:
        sif = ancient('LJA.sif'),
        reads = select_chry_reads
    output:
        assm = 'output/target_assembly/chry_reads/lja/{sample_info}_{sample}_{read_type}.{mapq}.k{kmer}-K{resolvek}/assembly.fasta',
    log:
        lja = 'log/output/target_assembly/chry_reads/{sample_info}_{sample}_{read_type}.{mapq}.k{kmer}-K{resolvek}.lja.log',
    benchmark:
        'rsrc/output/target_assembly/chry_reads/{sample_info}_{sample}_{read_type}.{mapq}.k{kmer}-K{resolvek}.lja.t24.rsrc',
    # conda:
    #     '../../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = CONSTRAINT_ALL_SAMPLES
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 131072 * attempt,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    params:
        out_folder = lambda wildcards, output: output.assm.rsplit('/', 1)[0]
    shell:
        'module load Singularity && singularity exec '
        '{input.sif} '
        'lja -o {params.out_folder} --reads {input.reads} -t {threads} -k {wildcards.kmer} -K {wildcards.resolvek} &> {log.lja} '
