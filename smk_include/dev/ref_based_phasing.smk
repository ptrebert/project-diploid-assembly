
include: '../environments.smk'

CHROMOSOMES = ['chr' + str(i) for i in range(1, 23)]

# must be placed under WORKING_DIRECTORY/references
REFERENCE_ASSEMBLY = 'GRCh38_HGSVC2_noalt'  # must have ".fasta" file name extension
#VARIANT_CALLS = 'HG002_GRCh38_1_22_v4.2_benchmark'

VARIANT_CALLS = config['variant_vcf']

#STRANDSEQ_PATH = '/home/local/work/data/lansdorp/benchmark_variants/run_folder/input/fastq/NA24385_bccrc_ilany-75pe_sseq'
#HIFI_READS_PATH = '/home/local/work/data/lansdorp/benchmark_variants/run_folder/input/fastq'

STRANDSEQ_PATH = config['strandseq_path']
HIFI_READS_PATH = config['hifi_path']


def collect_strandseq_libs(path):

    path = os.path.abspath(path)
    print('Collecting Strand-seq FASTQs from path: {}'.format(path))

    sseq_reads_id = os.path.split(path)[-1]
    libs = set()
    for fastq in os.listdir(path):
        if not fastq.endswith('.fastq.gz'):
            continue
        lib_id, mate_id = fastq.strip('.fastq.gz').rsplit('_', 1)
        if mate_id not in ['1', '2']:
            raise ValueError('Cannot parse SSEQ file name: {}'.format(fastq))
        libs.add(lib_id)
    print(sseq_reads_id, '({})'.format(len(libs)))
    return sseq_reads_id, sorted(libs)


def extract_hifi_reads_id(path):

    path = os.path.abspath(path)
    print('Collecting HiFi reads FASTQ from path: {}'.format(path))

    hifi_lib_id = None
    for fq in os.listdir(path):
        if not fq.endswith('.fastq.gz'):
            continue
        lib_id = fq.rsplit('.', 2)[0]
        if hifi_lib_id is not None:
            raise ValueError('Found more than one HiFi read input file: {} / {}'.format(hifi_lib_id, lib_id))
        hifi_lib_id = lib_id
    print(hifi_lib_id)
    return hifi_lib_id


def find_script_path(script_name, subfolder=''):
 
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


def load_fofn_file(input, prefix='', sep=' '):

    if not hasattr(input, 'fofn'):
        raise AttributeError('Input does not have "fofn" attribute: {}'.format(input))
    file_path = input.fofn
    if not os.path.isfile(file_path):
        file_list = ['FOFN-DRY-RUN']
    else:
        with open(file_path, 'r') as dump:
            file_list = sorted([l.strip() for l in dump.readlines()])
            assert file_list, 'Empty fofn file: {}'.format(file_path)
    file_list = prefix + sep.join(file_list)
    return file_list


SSEQ_READS_ID, SSEQ_READ_LIBS = collect_strandseq_libs(STRANDSEQ_PATH)

HIFI_READS_ID = extract_hifi_reads_id(HIFI_READS_PATH)

REF_CONFIG = {
    'reference': REFERENCE_ASSEMBLY,
    'variant_calls': VARIANT_CALLS,
    'sseq_reads': SSEQ_READS_ID,
    'sseq_libs': SSEQ_READ_LIBS,
    'hifi_reads': HIFI_READS_ID
}


PATH_REFERENCE_PHASING = '{reference}/{sseq_reads}/{hifi_reads}'

localrules: master_reference_phasing,
            write_breakpointr_config_file,
            write_strandphaser_config_file,
            write_phased_vcf_splits_fofn

rule master_reference_phasing:
    input:
        'output/integrative_phasing/{reference}/{sseq_reads}/{hifi_reads}/{variant_calls}.wh-phased.vcf.bgz'.format(**REF_CONFIG),
        'output/integrative_phasing/{reference}/{sseq_reads}/{hifi_reads}/{variant_calls}.wh-phased.vcf.bgz.tbi'.format(**REF_CONFIG),
        'output/statistics/phasing/{reference}/{sseq_reads}/{hifi_reads}/{variant_calls}.wh-phased.vcf.stats'.format(**REF_CONFIG)


rule generate_references_and_input:
    output:
        protected('references/{reference}.fasta'.format(**REF_CONFIG)),
        protected('references/{variant_calls}.vcf.gz'.format(**REF_CONFIG)),
        protected('input/{hifi_reads}/{hifi_reads}.fastq.gz'.format(**REF_CONFIG)),
        protected(
            expand(
                'input/{sseq_reads}/{sseq_library}_1.fastq.gz',
                sseq_reads=REF_CONFIG['sseq_reads'],
                sseq_library=REF_CONFIG['sseq_libs']
            )
        ),
        protected(
            expand(
                'input/{sseq_reads}/{sseq_library}_2.fastq.gz',
                sseq_reads=REF_CONFIG['sseq_reads'],
                sseq_library=REF_CONFIG['sseq_libs']
            )
        )


rule samtools_index_fasta:
    input:
        ancient('{filepath}.fasta')
    output:
        '{filepath}.fasta.fai'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
        'samtools faidx {input}'


rule bgzip_file_copy:
    input:
        '{filepath}'
    output:
        '{filepath}.bgz'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
        'bgzip -c {input} > {output}'


rule bcftools_index_bgzipped_file:
    input:
        '{filepath}.bgz'
    output:
        '{filepath}.bgz.tbi'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
        'bcftools index --tbi {input}'


rule create_assembly_sequence_files:
    input:
        'references/{reference}.fasta.fai'
    output:
        directory('references/{reference}/sequences')
    run:
        import os
        output_dir = output[0]
        os.makedirs(output_dir, exist_ok=True)
        with open(input[0], 'r') as fai:
            for line in fai:
                seq_name = line.split('\t')[0]
                if seq_name not in CHROMOSOMES:
                    continue
                output_path = os.path.join(output_dir, seq_name + '.seq')
                if os.path.isfile(output_path):
                    with open(output_path, 'r') as seq_file:
                        seq_entry = seq_file.read().strip()
                        if seq_entry == line.strip():
                            continue
                with open(output_path, 'w') as dump:
                    _ = dump.write(line)


rule generate_bwa_index:
    input:
        reference = ancient('references/{reference}.fasta')
    output:
        multiext(
            'references/{reference}/bwa_index/{reference}',
            '.amb', '.ann', '.bwt', '.pac', '.sa'
        )
    log:
        'log/references/bwa_index/{reference}.log'
    benchmark:
        'run/references/bwa_index/{reference}.rsrc'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    conda:
         '../../environment/conda/conda_biotools.yml'
    params:
        prefix = lambda wildcards, output: output[0].rsplit('.', 1)[0]
    shell:
        'bwa index -p {params.prefix} {input.reference} &> {log}'


rule pbmm2_reads_to_reference_alignment_pacbio_fastq:
    input:
        reads = ancient('input/{hifi_reads}/{hifi_reads}.fastq.gz'),
        reference = ancient('references/{reference}.fasta')
    output:
        bam = 'output/alignments/reads_to_reference/{hifi_reads}_map-to_{reference}.psort.sam.bam',
    log:
        'log/output/alignments/reads_to_reference/{hifi_reads}_map-to_{reference}.pbmm2.log'
    benchmark:
        'run/output/alignments/reads_to_reference/{{hifi_reads}}_map-to_{{reference}}.pbmm2.t{}.rsrc'.format(config['num_cpu_high'])
    conda:
        '../../environment/conda/conda_pbtools.yml'
    wildcard_constraints:
        hifi_reads = HIFI_READS_ID,
        reference = REF_CONFIG['reference']
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((110592 if attempt <= 1 else 188416) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 110592 if attempt <= 1 else 188416,
        runtime_hrs = lambda wildcards, attempt: 8 * attempt
    params:
        align_threads = config['num_cpu_high'] - 2,
        sort_threads = 2,
        preset = ' CCS --min-length 5000 ',  # hard-code CCS preset
        individual = lambda wildcards: wildcards.hifi_reads.split('_')[0],
        tempdir = lambda wildcards: os.path.join(
                                        'temp', 'pbmm2', 
                                        wildcards.hifi_reads,
                                        wildcards.reference
                                    )
    shell:
        'TMPDIR={params.tempdir} '
        'pbmm2 align --log-level INFO --sort --sort-memory {resources.mem_per_cpu_mb}M --no-bai '
            ' --alignment-threads {params.align_threads} --sort-threads {params.sort_threads} '
            ' --preset {params.preset} --rg "@RG\\tID:1\\tSM:{params.individual}" --sample {params.individual} '
            ' {input.reference} {input.reads} {output.bam} &> {log}'


rule samtools_index_bam_alignment:
    """
    - multi-threaded index generation seems to swallow error (e.g., bad blocks / error 33)
    - WH claimed that, by experience, single-threaded is often faster
    """
    input:
        bam = '{filepath}.bam'
    output:
        bai = '{filepath}.bam.bai'
    benchmark:
        'run/{filepath}.idx-bai.t1.rsrc'
    resources:
        runtime_hrs = lambda wildcards, attempt: 0 if attempt <= 1 else 16 * attempt,
        runtime_min = lambda wildcards, attempt: 20  # short time limit for Strand-seq alignments
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
        "samtools index {input.bam}"


rule bwa_strandseq_to_reference_alignment:
    input:
        mate1 = ancient('input/{sseq_reads}/{individual}_{sample_id}_1.fastq.gz'),
        mate2 = ancient('input/{sseq_reads}/{individual}_{sample_id}_2.fastq.gz'),
        ref_index = 'references/{reference}/bwa_index/{reference}.bwt',
    output:
        bam = 'output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/aln/{individual}_{sample_id}.filt.sam.bam'
    log:
        bwa = 'log/output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/{individual}_{sample_id}.bwa.log',
        samtools = 'log/output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/{individual}_{sample_id}.samtools.log',
    benchmark:
        os.path.join('run/output/alignments/strandseq_to_reference/{reference}/{sseq_reads}',
                     '{individual}_{sample_id}' + '.t{}.rsrc'.format(config['num_cpu_low']))
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(12288 * attempt / config['num_cpu_low']),
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt
    params:
        idx_prefix = lambda wildcards, input: input.ref_index.rsplit('.', 1)[0],
        discard_flag = config['bwa_strandseq_aln_discard']
    shell:
        'bwa mem -t {threads}'
            ' -R "@RG\\tID:{wildcards.individual}_{wildcards.sample_id}\\tPL:Illumina\\tSM:{wildcards.individual}"'
            ' -v 2 {params.idx_prefix} {input.mate1} {input.mate2} 2> {log.bwa} | '
            ' samtools view -b -F {params.discard_flag} /dev/stdin > {output.bam} 2> {log.samtools}'


rule samtools_position_sort_strandseq_reads:
    """
    Since Strand-seq alignments are small, make dedicated
    samtools sort rule with lower resource requirements
    """
    input:
        'output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/aln/{sseq_library}.filt.sam.bam'
    output:
        temp('output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/sort/{sseq_library}.filt.psort.sam.bam')
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sseq_library = '[A-Za-z0-9\-_]+'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_total_mb = lambda wildcards, attempt: config['num_cpu_low'] * 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'samtools sort -m {resources.mem_per_cpu_mb}M --threads {threads} -o {output} {input}'


rule mark_duplicate_reads_strandseq:
    input:
        'output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/sort/{sseq_library}.filt.psort.sam.bam'
    output:
        'output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/{sseq_library}.filt.psort.mdup.sam.bam'
    log:
        'log/output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/{sseq_library}.filt.psort.mdup.log'
    benchmark:
        'run/output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/{sseq_library}.filt.psort.mdup.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sseq_library = '[A-Za-z0-9\-_]+'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_total_mb = lambda wildcards, attempt: config['num_cpu_low'] * 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'sambamba markdup -t {threads} --overflow-list-size 600000 {input} {output} &> {log}'


rule write_breakpointr_config_file:
    """
    As long as aggregate-style input functions downstream of
    checkpoints are problematic in a cluster environment,
    avoid complications by making the config writing local
    (aggregate should work), and explicitly write a file
    containing just the input folder for breakpointR
    """
    input:
        setup_ok = rules.install_rlib_breakpointr.output.check,
        reference = ancient('references/{reference}.fasta'),
        bam = expand(
            'output/alignments/strandseq_to_reference/{{reference}}/{{sseq_reads}}/{sseq_library}.filt.psort.mdup.sam.bam',
            sseq_library=REF_CONFIG['sseq_libs']
        ),
        bai = expand(
            'output/alignments/strandseq_to_reference/{{reference}}/{{sseq_reads}}/{sseq_library}.filt.psort.mdup.sam.bam.bai',
            sseq_library=REF_CONFIG['sseq_libs']
        )
    output:
        cfg = 'output/integrative_phasing/processing/config_files/{reference}/{sseq_reads}/breakpointr.config',
        input_dir = 'output/integrative_phasing/processing/config_files/{reference}/{sseq_reads}/breakpointr.input',
    threads: 1
    params:
        bp_cpu = config['num_cpu_high']
    run:
        outfolder = os.path.dirname(input.bam[0])

        config_rows = [
            '[General]',
            'numCPU = ' + str(params.bp_cpu),  # due to a bug in breakpointr, this value has to be repeated on the CLI
            'reuse.existing.files = FALSE',
            '',
            '[breakpointR]',
            'windowsize = 500000',
            'binMethod = "size"',
            'pairedEndReads = TRUE',
            'pair2frgm = FALSE',
            'min.mapq = 10',
            'filtAlt = TRUE',
            'background = 0.1',
            'minReads = 50'
        ]

        with open(output.cfg, 'w') as dump:
            _ = dump.write('\n'.join(config_rows) + '\n')

        with open(output.input_dir, 'w') as dump:
            _ = dump.write(outfolder + '\n')


rule run_breakpointr:
    input:
        cfg = rules.write_breakpointr_config_file.output.cfg,
        fofn = rules.write_breakpointr_config_file.output.input_dir
    output:
        wc_reg = 'output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/{reference}.WCregions.txt',
        cfg = 'output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/run/breakpointR.config',
        rdme = 'output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/run/README.txt',
        bps = directory('output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/run/breakpoints'),
        bwf = directory('output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/run/browserfiles'),
        data = directory('output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/run/data'),
        plots = directory('output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/run/plots')
    log:
        'log/output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/breakpointr.log'
    benchmark:
        os.path.join('run/output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}',
                     'breakpointr.t{}.rsrc'.format(config['num_cpu_high'])
                     )
    conda:
        '../../environment/conda/conda_rscript.yml'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(49152    * attempt / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 12 * attempt
    params:
        output_dir = lambda wildcards, output: os.path.dirname(output.rdme),
        input_dir = lambda wildcards, input: load_fofn_file(input),
        script_exec = lambda wildcards: find_script_path('run_breakpointr.R')
    shell:
        '{params.script_exec} {params.input_dir} {input.cfg} {params.output_dir} {threads} {output.wc_reg} &> {log}'


rule unphase_reference_vcf:
    input:
        vcf = ancient('references/{variant_calls}.vcf.gz'),
    output:
        vcf = 'references/{variant_calls}.unphased.vcf'
    log:
        'log/references/{variant_calls}.unphased.log'
    benchmark:
        'run/references/{variant_calls}.unphased.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        variant_calls = REF_CONFIG['variant_calls']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'whatshap --debug unphase {input.vcf} > {output.vcf} 2> {log}'


rule write_strandphaser_config_file:
    """
    As long as aggregate-style input functions downstream of
    checkpoints are problematic in a cluster environment,
    avoid complications by making the config writing local
    (aggregate should work), and explicitly write a file
    containing just the input folder for StrandPhaseR
    """
    input:
        setup_ok = rules.install_rlib_strandphaser.output.check,
        reference = ancient('references/{reference}.fasta'),
        bam = expand(
            'output/alignments/strandseq_to_reference/{{reference}}/{{sseq_reads}}/{sseq_library}.filt.psort.mdup.sam.bam',
            sseq_library=REF_CONFIG['sseq_libs']
        ),
        bai = expand(
            'output/alignments/strandseq_to_reference/{{reference}}/{{sseq_reads}}/{sseq_library}.filt.psort.mdup.sam.bam.bai',
            sseq_library=REF_CONFIG['sseq_libs']
        )
    output:
        cfg = 'output/integrative_phasing/processing/config_files/{reference}/{sseq_reads}/strandphaser.config',
        input_dir = 'output/integrative_phasing/processing/config_files/{reference}/{sseq_reads}/strandphaser.input',
    threads: 1
    params:
        sp_cpu = config['num_cpu_high']
    run:
        import os

        outfolder = os.path.dirname(input.bam[0])

        config_rows = [
            '[General]',
            'numCPU = ' + str(params.sp_cpu),
            'pairedEndReads = TRUE',
            'min.mapq = 10',
            '',
            '[StrandPhaseR]',
            'min.baseq = 20',
            'num.iterations = 2',
            'translateBases = TRUE',
            'splitPhasedReads = TRUE'
        ]

        with open(output.cfg, 'w') as dump:
            _ = dump.write('\n'.join(config_rows) + '\n')

        with open(output.input_dir, 'w') as dump:
            _ = dump.write(outfolder + '\n')


rule run_strandphaser:
    input:
        cfg = 'output/integrative_phasing/processing/config_files/{reference}/{sseq_reads}/strandphaser.config',
        fofn = 'output/integrative_phasing/processing/config_files/{reference}/{sseq_reads}/strandphaser.input',
        wc_regions = 'output/integrative_phasing/processing/breakpointr/{reference}/{sseq_reads}/{reference}.WCregions.txt',
        variant_calls = 'references/{variant_calls}.unphased.vcf'
    output:
        browser = directory('output/integrative_phasing/processing/strandphaser/' + PATH_REFERENCE_PHASING + '/{variant_calls}/browserFiles'),
        data = directory('output/integrative_phasing/processing/strandphaser/' + PATH_REFERENCE_PHASING + '/{variant_calls}/data'),
        phased = directory('output/integrative_phasing/processing/strandphaser/' + PATH_REFERENCE_PHASING + '/{variant_calls}/Phased'),
        maps = directory('output/integrative_phasing/processing/strandphaser/' + PATH_REFERENCE_PHASING + '/{variant_calls}/SingleCellHaps'),
        vcf_dir = directory('output/integrative_phasing/processing/strandphaser/' + PATH_REFERENCE_PHASING + '/{variant_calls}/VCFfiles'),
        cfg = 'output/integrative_phasing/processing/strandphaser/' + PATH_REFERENCE_PHASING + '/{variant_calls}/StrandPhaseR.config',
    log:
        stp = 'log/output/integrative_phasing/processing/strandphaser/' + PATH_REFERENCE_PHASING + '/{variant_calls}.phased.log',
    benchmark:
        os.path.join('run/output/integrative_phasing/processing/strandphaser',
                     PATH_REFERENCE_PHASING + '/{variant_calls}' + '.phased.t{}.rsrc'.format(config['num_cpu_high'])
                     )
    conda:
        '../../environment/conda/conda_rscript.yml'
    wildcard_constraints:
        sample = HIFI_READS_ID,
        reference = REF_CONFIG['reference'],
        variant_calls = REF_CONFIG['variant_calls']
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(49152 * attempt / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 49152 * attempt,
        runtime_hrs = lambda wildcards, attempt: 12 * attempt
    params:
        input_dir = lambda wildcards, input: load_fofn_file(input),
        output_dir = lambda wildcards, output: os.path.dirname(output.cfg),
        individual = lambda wildcards: wildcards.sseq_reads.split('_')[0],
        script_exec = lambda wildcards: find_script_path('run_strandphaser.R')
    shell:
        '{params.script_exec} {params.input_dir} {input.cfg} '
            ' {input.variant_calls} {input.wc_regions} '
            ' {params.output_dir} {params.individual} &> {log.stp}'


rule write_strandphaser_split_vcf_fofn:
    input:
        rules.run_strandphaser.output.vcf_dir
    output:
        fofn = 'output/integrative_phasing/processing/strandphaser/' + PATH_REFERENCE_PHASING + '/{variant_calls}.spr-phased.fofn'
    resources:
        mem_total_mb = 2048,
        mem_per_cpu_mb = 2048,
    run:
        # Note that this rule does not call an aggregate function because
        # run_strandphaser is not a checkpoint... I *think* the problem
        # was the downstream call "run_integrative_phasing", where direct
        # access to an individual VCF in /VCFfiles/{sequence}.vcf did not work
        # See comment below for rule "run_integrative_phasing"

        import os

        input_dir = input[0]
        input_vcfs = sorted([f for f in os.listdir(input_dir) if f.endswith('.vcf')])
        if len(input_vcfs) == 0:
            raise RuntimeError('No StrandPhaseR-phased VCFs in /VCFfiles output dir. '
                               'Cannot create fofn file: {}'.format(output.fofn))

        with open(output.fofn, 'w') as fofn:
            for file_name in input_vcfs:
                file_path = os.path.join(input_dir, file_name)
                if not os.path.isfile(file_path):
                    if os.path.isdir(file_path):
                        # definitely wrong...
                        raise RuntimeError('Found directory, but expected VCF file: {} '
                                           '(write fofn file: {}'.format(file_path, output.fofn))
                _ = fofn.write(file_path + '\n')


rule merge_strandphaser_output:
    input:
        fofn = 'output/integrative_phasing/processing/strandphaser/' + PATH_REFERENCE_PHASING + '/{variant_calls}.spr-phased.fofn'
    output:
        vcf = 'output/integrative_phasing/' + PATH_REFERENCE_PHASING + '/{variant_calls}.spr-phased.vcf'
    log:
        'log/output/integrative_phasing/' + PATH_REFERENCE_PHASING + '/{variant_calls}.concat-phased.log'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bcftools concat -f {input.fofn} --output-type v --output {output} &> {log}'


rule run_integrative_phasing:
    """
    Why this complicated command line?
    One of those "why-is-this-a-problem" situations; for some reason, Snakemake
    does not recognize the output of the rule "run_strandphaser"
    as the producer for the individual VCF splits. So, merge, and split again...
    """
    input:
        vcf = 'references/{variant_calls}.unphased.vcf.bgz',
        tbi = 'references/{variant_calls}.unphased.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/{hifi_reads}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/{hifi_reads}_map-to_{reference}.psort.sam.bam.bai',
        fasta = ancient('references/{reference}.fasta'),
        fai = 'references/{reference}.fasta.fai',
        spr_phased = 'output/integrative_phasing/' + PATH_REFERENCE_PHASING + '/{variant_calls}.spr-phased.vcf'
    output:
        vcf = 'output/integrative_phasing/processing/whatshap/' + PATH_REFERENCE_PHASING + '/{variant_calls}.{sequence}.phased.vcf'
    log:
        'log/output/integrative_phasing/processing/whatshap/' + PATH_REFERENCE_PHASING + '/{variant_calls}.{sequence}.phased.log'
    benchmark:
        'run/output/integrative_phasing/processing/whatshap/' + PATH_REFERENCE_PHASING + '/{variant_calls}.{sequence}.phased.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sample = HIFI_READS_ID,
        reference = REF_CONFIG['reference'],
        variant_calls = REF_CONFIG['variant_calls']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'whatshap --debug phase --chromosome {wildcards.sequence} --ignore-read-groups --indels '
            '--reference {input.fasta} {input.vcf} {input.bam} {input.spr_phased} 2> {log} '
            ' | egrep "^(#|{wildcards.sequence}\s)" > {output}'


rule write_phased_vcf_splits_fofn:
    input:
        vcf_splits = expand(
            'output/integrative_phasing/processing/whatshap/{{reference}}/{{sseq_reads}}/{{hifi_reads}}/{{variant_calls}}.{sequence}.phased.vcf',
            sequence=CHROMOSOMES
        )
    output:
        fofn = 'output/integrative_phasing/processing/whatshap/' + PATH_REFERENCE_PHASING + '/{variant_calls}.wh-phased.fofn'
    run:
        with open(output.fofn, 'w') as dump:
            for file_path in sorted(input.vcf_splits):
                if not os.path.isfile(file_path):
                    if os.path.isdir(file_path):
                        # this is definitely wrong
                        raise AssertionError('Expected file path for INTPHASE VCF split merge, but received directory: {}'.format(file_path))
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule strandseq_dga_merge_sequence_phased_vcf_files:
    input:
        fofn = rules.write_phased_vcf_splits_fofn.output.fofn
    output:
        'output/integrative_phasing/' + PATH_REFERENCE_PHASING + '/{variant_calls}.wh-phased.vcf'
    log:
        'log/output/integrative_phasing/' + PATH_REFERENCE_PHASING + '/{variant_calls}.wh-phased.concat.log'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
             mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
             mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    shell:
        'bcftools concat -f {input.fofn} --output {output} --output-type v &> {log}'


rule compute_whatshap_phased_vcf_stats:
    input:
        vcf = 'output/integrative_phasing/' + PATH_REFERENCE_PHASING + '/{variant_calls}.wh-phased.vcf.bgz',
        idx = 'output/integrative_phasing/' + PATH_REFERENCE_PHASING + '/{variant_calls}.wh-phased.vcf.bgz.tbi'
    output:
        bcf_stats = 'output/statistics/phasing/' + PATH_REFERENCE_PHASING + '/{variant_calls}.wh-phased.vcf.stats',
        wh_stats_tsv = 'output/statistics/phasing/' + PATH_REFERENCE_PHASING + '/{variant_calls}.wh-phased.stats.tsv',
        wh_stats_txt = 'output/statistics/phasing/' + PATH_REFERENCE_PHASING + '/{variant_calls}.wh-phased.stats.txt'
    log:
        'log/output/statistics/phasing/' + PATH_REFERENCE_PHASING + '/{variant_calls}.wh-phased.stats.log'
    conda:
        '../../environment/conda/conda_biotools.yml'
    priority: 200
    resources:
             mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
             mem_total_mb = lambda wildcards, attempt: 4096 * attempt
    shell:
        'bcftools stats {input.vcf} > {output.bcf_stats} 2> {log}'
            ' && '
            'whatshap stats --tsv {output.wh_stats_tsv} {input.vcf} > {output.wh_stats_txt} 2>> {log}'


