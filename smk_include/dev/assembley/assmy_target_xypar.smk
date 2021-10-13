

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
    import pathlib as pl
    filename = '{sample_info}_{sample}_{read_type}.chrY-reads.{mapq}.'.format(**wildcards)
    reads = list(pl.Path('output/read_subsets/chry').glob(f'{filename}*.gz'))
    assert len(reads) == 1
    return reads[0]


rule hifiasm_chry_targeted_assembly:
    input:
        reads = select_chry_reads
    output:
        primary_contigs = 'output/target_assembly/chry_reads/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.p_ctg.gfa',
        raw_unitigs = 'output/target_assembly/chry_reads/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.r_utg.gfa',
    log:
        hifiasm = 'log/output/target_assembly/chry_reads/{sample_info}_{sample}_{read_type}.{mapq}.hifiasm.log',
    benchmark:
        'rsrc/output/target_assembly/chry_reads/{sample_info}_{sample}_{read_type}.{mapq}.hifiasm.t12.rsrc',
    conda:
        '../../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = CONSTRAINT_ALL_SAMPLES
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        prefix = lambda wildcards, output: output.primary_contigs.rsplit('.', 2)[0],
    shell:
        'hifiasm -o {params.prefix} -t {threads} --write-ec --write-paf --primary {input.reads} &> {log.hifiasm}'


rule mbg_chry_targeted_assembly:
    """
    sif = ancient('mbg.master.sif'),
    """
    input:
        sif = ancient('mbg.UnitigResolve.sif'),
        reads = select_chry_reads
    output:
        graph = 'output/target_assembly/chry_reads/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.MBG-k{kmer}-w{window}.gfa',
        paths = 'output/target_assembly/chry_reads/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.MBG-k{kmer}-w{window}.gaf',
    log:
        'log/output/target_assembly/chry_reads/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.MBG-k{kmer}-w{window}.log',
    benchmark:
        'rsrc/output/target_assembly/chry_reads/mbg/{sample}/{sample_info}_{sample}_{read_type}.{mapq}.MBG-k{kmer}-w{window}.rsrc',
    wildcard_constraints:
        read_type = 'HIFIEC',
        sample = '(NA193N7|NA193NN|AFR4MIX)'
#    conda:
#        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 65536 + 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    shell:
        'module load Singularity && singularity exec '
        '--bind /:/hilbert {input.sif} '
        'MBG -i {input.reads} -t {threads} '
            '-k {wildcards.kmer} -w {wildcards.window} '
            '--error-masking no --include-end-kmers '
            '--out {output.graph} --output-sequence-paths {output.paths} &> {log}'
