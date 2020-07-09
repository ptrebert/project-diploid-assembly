
include: 'prep_custom_references.smk'

QVEST_CONFIG = {
    'discard_flag': 2816,  # not primary OR qc fail OR supplementary
    'cov_threshold': '{mean} + 3 {stddev}',  # avoid that freebayes gets tangled in ultra high cov regions
    'genome_size': int(3.1e9),  # this is roughly the sequence length of the HGSVC2 reference
    'ref_assembly': 'GRCh38_HGSVC2_noalt',
    'skip_short_read_sources': ['PRJEB3381', 'PRJEB9396']
}


def find_sample_short_reads(sample):

    data_sources = config['sample_description_' + sample]['data_sources']
    short_reads = []
    for readset_record in data_sources:
        for readset_type, readset_desc in readset_record.items():
            if readset_type != 'short_reads':
                continue
            readset_sample, readset_name = readset_desc['readset'].split('_', 1)
            assert readset_sample == sample, 'Sample mismatch: {} / {}'.format(sample, readset_desc)

            bioproject = None
            try:
                bioproject = readset_desc['bioproject']
            except KeyError:
                pass
            if bioproject in KMER_CONFIG['skip_short_read_sources']:
                continue

            short_reads.append(readset_name)

    if not short_reads:
        raise ValueError('No short read data available for sample {}'.format(sample))
    return short_reads


def determine_possible_computations(wildcards):
    """
    NA19239_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta
    """
    module_outputs = {
        'bam_stats': os.path.join('output/alignments/short_to_phased_assembly',
                     '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.mdup.stats')
    }

    fix_wildcards = {
        'known_ref': QVEST_CONFIG['ref_assembly']
    }

    compute_results = set()

    search_path = os.path.join(os.getcwd(), 'output/evaluation/phased_assemblies')
    for ps_assm in os.listdir(search_path):
        if not ps_assm.endswith('.fasta'):
            continue
        assm_base, hap, polisher, ext = ps_assm.split('.')
        sample, assm_reads = assm_base.split('_', 1)
        tmp = dict(fix_wildcards)
        tmp['sample'] = sample
        tmp['polisher'] = polisher
        tmp['assembly'] = assm_reads
        tmp['hap'] = hap
        short_reads = find_sample_short_reads(sample)
        for sr in short_reads:
            tmp['readset'] = sr
            for target in module_outputs.values():
                fmt_target = target.format(**tmp)
                compute_results.add(fmt_target)
    
    return sorted(compute_results)


localrules: master_qv_estimate


rule master_qv_estimate:
    input:
        determine_possible_computations


rule bwa_short_to_haploid_assembly_alignment:
    input:
        mate1 = 'input/fastq/{individual}_{library_id}_short_1.fastq.gz',
        mate2 = 'input/fastq/{individual}_{library_id}_short_2.fastq.gz',
        ref_index = os.path.join(
            'output/evaluation/phased_assemblies',
            '{individual}_{assembly}.{hap}.{polisher}',
            'bwa_index',
            '{individual}_{assembly}.{hap}.{polisher}.bwt'
        )
    output:
        bam = os.path.join(
            'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.filt.sort.sam.bam'
        )
    log:
        bwa = os.path.join(
            'log/output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.filt.sort.bwa.log'),
        st_view = os.path.join(
            'log/output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.filt.sort.st-view.log'),
        st_sort = os.path.join(
            'log/output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.filt.sort.st-sort.log'),
    benchmark:
        os.path.join(
            'run/output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.filt.sort' + '.t{}.rsrc'.format(config['num_cpu_high']))
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((24768 + 12288 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 24768 + 12288 * attempt,
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    params:
        idx_prefix = lambda wildcards, input: input.ref_index.rsplit('.', 1)[0],
        discard_flag = QVEST_CONFIG['discard_flag'],
        sort_mem_mb = 4096,
        sort_threads = 6
    shell:
        'bwa mem -t {threads}'
            ' -R "@RG\\tID:{wildcards.individual}_{wildcards.library_id}\\tPL:Illumina\\tSM:{wildcards.individual}"'
            ' -v 2 {params.idx_prefix} {input.mate1} {input.mate2} 2> {log.bwa} | '
            ' samtools view -b -F {params.discard_flag} /dev/stdin 2> {log.st_view} | '
            ' samtools sort -l 6 -@ {params.sort_threads} -m {params.sort_mem_mb}M -O BAM -o {output} /dev/stdin 2> {log.st_sort}'


rule mark_duplicate_reads:
    input:
        bam = os.path.join(
            'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.filt.sort.sam.bam'
        )
    output:
        bam = os.path.join(
            'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup.sam.bam'
        )
    log:
        os.path.join(
            'log', 'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup.log'
        )
    benchmark:
        os.path.join(
            'run', 'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup' + '.t{}.rsrc'.format(config['num_cpu_low'])
        )
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(12288 * attempt / config['num_cpu_low']),
        mem_total_mb = lambda wildcards, attempt: 12288 * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    shell:
        'sambamba markdup -t {threads} '
        '--sort-buffer-size=8092 '
        '--overflow-list-size 1000000 '  # default 200 000 ; increase to avoid "too many open files" issue
        '--hash-table-size 524288 '  # default 262 144 ; increase to avoid "too many open files" issue
        '{input} {output} 2> {log}'


rule compute_alignments_stats:
    input:
        bam = os.path.join(
            'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup.sam.bam'
        ),
        bai = os.path.join(
            'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup.sam.bam.bai'
        )
    output:
        os.path.join(
            'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup.stats'
        )
    log:
        os.path.join(
            'log', 'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup.stats.log'
        )
    benchmark:
        os.path.join(
            'run', 'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup.stats' + '.t{}.rsrc'.format(config['num_cpu_low'])
        )
    conda:
         '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int(4096 * attempt / config['num_cpu_low']),
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    shell:
        'samtools stats --remove-dups --threads {threads} {input.bam} > {output} 2> {log}'
