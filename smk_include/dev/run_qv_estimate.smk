
include: 'prep_custom_references.smk'

QVEST_CONFIG = {
    'discard_flag': 2816,  # not primary OR qc fail OR supplementary
    'cov_threshold': '{mean} + 3 * {stddev}',  # avoid that freebayes gets tangled in ultra high cov regions
    'genome_size': int(3.1e9),  # this is roughly the sequence length of the HGSVC2 reference
    'freebayes_timeout_sec': 14400,
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
            if bioproject in QVEST_CONFIG['skip_short_read_sources']:
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
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.mdup.stats'),
        'raw_calls': os.path.join('output/evaluation/qv_estimation/variant_stats/00-raw',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.stats'),
        'hom_snps': os.path.join('output/evaluation/qv_estimation/variant_stats/30-split-gtype',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.snps.hom.stats'),
        'hom_indels': os.path.join('output/evaluation/qv_estimation/variant_stats/30-split-gtype',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.indels.hom.stats'),
        'het_snps': os.path.join('output/evaluation/qv_estimation/variant_stats/30-split-gtype',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.snps.het.stats'),
        'het_indels': os.path.join('output/evaluation/qv_estimation/variant_stats/30-split-gtype',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.indels.het.stats'),
    }

    fix_wildcards = {
        'known_ref': QVEST_CONFIG['ref_assembly']
    }

    compute_results = set()

    search_path = os.path.join(os.getcwd(), 'output/evaluation/phased_assemblies')
    for ps_assm in os.listdir(search_path):
        if ps_assm.startswith('v1'):
            version, new_name = ps_assm.split('_', 1)
            os.rename(os.path.join(search_path, ps_assm), os.path.join(search_path, new_name))
            assm_file = new_name
        else:
            assm_file = ps_assm
        if not assm_file.endswith('.fasta'):
            continue
        assm_base, hap, polisher, ext = assm_file.split('.')
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


localrules: master_qv_estimate,
            compute_coverage_statistics,
            convert_faidx_to_regions


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
    wildcard_constraints:
        individual = '[A-Z0-9]+'
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


def weighted_avg_and_std(values, weights):
    """
    https://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
    Return the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    """
    import numpy as np
    try:
        average = np.average(values, weights=weights)
    except ZeroDivisionError:
        return 0, 0
    variance = np.average((values-average)**2, weights=weights)
    return average, np.sqrt(variance)


rule compute_coverage_statistics:
    input:
        stats_file = os.path.join('output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup.stats'),
        faidx = os.path.join('output/evaluation/phased_assemblies',
            '{individual}_{assembly}.{hap}.{polisher}.fasta.fai')
    output:
        os.path.join('output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup.cov.tsv'),
    wildcard_constraints:
        individual = '[A-Z0-9]+'
    run:
        import numpy as np

        with open(input.faidx, 'r') as faidx:
            seqlen = sum([int(line.split()[1]) for line in faidx.readlines()])

        # coverage range from [ 0 ... 1001 ]
        coverages = np.zeros(1002, dtype=np.int32)

        with open(input.stats_file, 'r') as stats:
            for line in stats:
                if not line.startswith('COV'):
                    continue
                _, cov_range, cov, count = line.strip().split()
                if '<' in cov_range:
                    cov = 1001
                else:
                    cov = int(cov)
                coverages[cov] += int(count)
        
        # number of bases with zero coverage
        coverages[0] = seqlen - coverages.sum()

        # based on some initial tests, seems to be unlikely,
        # but in principle many coverage values could be zero
        observed_coverage_values = np.array(coverages > 0, dtype=np.bool)
        cov_values = np.array(range(1002), dtype=np.int32)[observed_coverage_values]
        cov_counts = coverages[observed_coverage_values]
        avg, std = weighted_avg_and_std(cov_values, cov_counts)

        with open(output[0], 'w') as dump:
            _ = dump.write('cov_mean\t{}\n'.format(round(avg, 2)))
            _ = dump.write('cov_stddev\t{}\n'.format(round(std, 2)))
            _ = dump.write('cov_bp_nonzero\t{}\n'.format(coverages[1:].sum()))
            _ = dump.write('cov_bp_zero\t{}\n'.format(coverages[0]))
            _ = dump.write('seq_length\t{}\n'.format(seqlen))


rule convert_faidx_to_regions:
    """
    Regions file needed to run freebayes parallel script
    """
    input:
        faidx = os.path.join('output/evaluation/phased_assemblies',
            '{individual}_{assembly}.{hap}.{polisher}.fasta.fai')
    output:
        regions = os.path.join('output/evaluation/phased_assemblies',
            '{individual}_{assembly}.{hap}.{polisher}.fasta.regions')
    wildcard_constraints:
        individual = '[A-Z0-9]+'
    run:
        regions = []
        with open(input.faidx, 'r') as faidx:
            for line in faidx:
                parts = line.split()
                seq_id, seq_size = parts[:2]
                regions.append('{}:0-{}'.format(seq_id, seq_size))
        
        with open(output.regions, 'w') as dump:
            _ = dump.write('\n'.join(regions) + '\n')


def compute_coverage_limit(cov_file):

    # coverages for Illumina data varies wildly, and freebayes' runtime
    # seems to be very sensitive to coverage. Set some
    # hard limit no matter what the observed coverage actually is

    HARD_LIMIT = 200

    if not os.path.isfile(cov_file):
        # assume dry run
        return HARD_LIMIT

    with open(cov_file, 'r') as cov_info:
        _, cov_mean = cov_info.readline().split()
        _, cov_stddev = cov_info.readline().split()

    cov_limit = eval(QVEST_CONFIG['cov_threshold'].format(**{'mean': cov_mean, 'stddev': cov_stddev}))
    cov_limit = int(round(cov_limit, 0))
    cov_limit = max(min(cov_limit, HARD_LIMIT), 1)  # max 1: not sure if freebayes would handle 0 correctly
    return cov_limit


rule qvest_freebayes_call_variants:
    input:
        bam = os.path.join(
            'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup.sam.bam'),
        bai = os.path.join(
            'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup.sam.bam.bai'),
        reference = os.path.join(
            'output/evaluation/phased_assemblies',
            '{individual}_{assembly}.{hap}.{polisher}.fasta'),
        ref_idx = os.path.join(
            'output/evaluation/phased_assemblies',
            '{individual}_{assembly}.{hap}.{polisher}.fasta.fai'),
        ref_regions = os.path.join(
            'output/evaluation/phased_assemblies',
            '{individual}_{assembly}.{hap}.{polisher}.fasta.regions'),
        cov_file = os.path.join(
            'output/alignments/short_to_phased_assembly',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.mdup.cov.tsv')
    output:
        os.path.join(
            'output/evaluation/qv_estimation/variant_calls/00-raw',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.vcf.bgz'
            )
    log:
        os.path.join(
            'log', 'output/evaluation/qv_estimation/variant_calls/00-raw',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.fb-call.log'
            )
    benchmark:
        os.path.join(
            'run', 'output/evaluation/qv_estimation/variant_calls/00-raw',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.fb-call' + '.t{}.rsrc'.format(config['num_cpu_high'])
            )
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        individual = '[A-Z0-9]+'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((32768 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    params:
        timeout = QVEST_CONFIG['freebayes_timeout_sec'],
        script_exec = lambda wildcards: find_script_path('fb-parallel-timeout.sh'),
        skip_cov = lambda wildcards, input: compute_coverage_limit(input.cov_file)
    shell:
        '{params.script_exec} {input.ref_regions} {threads} {params.timeout} {log}'
            ' --use-best-n-alleles 4 --strict-vcf --fasta-reference {input.reference}'
            ' --skip-coverage {params.skip_cov} {input.bam} | bgzip -c /dev/stdin > {output}'


rule variant_calls_qfilter_biallelic:
    input:
        vcf = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/00-raw',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.vcf.bgz'
            ),
        tbi = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/00-raw',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.vcf.bgz.tbi'
            ),
        ref = os.path.join(
            'output/evaluation/phased_assemblies',
            '{individual}_{assembly}.{hap}.{polisher}.fasta'),
    output:
        vcf = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/10-qfilter-biallelic',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.q10-bia.vcf.bgz'
            ),
    wildcard_constraints:
        individual = '[A-Z0-9]+'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bcftools filter --output-type u --include "QUAL>=10" | '
        ' --min-alleles 2 --max-alleles 2 {input.vcf} | '
        ' bcftools norm --fasta-ref {input.ref} --output /dev/stdout --output-type v | '
        ' bgzip -c /dev/stdin > {output}'


rule split_callset_by_variant_type:
    input:
        vcf = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/10-qfilter-biallelic',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.q10-bia.vcf.bgz'
            ),
        tbi = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/10-qfilter-biallelic',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.q10-bia.vcf.bgz.tbi'
            ),
    output:
        snps = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/20-split-vtype',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.snps.vcf.bgz'
            ),
        indels = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/20-split-vtype',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.indels.vcf.bgz'
            ),
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
         'bcftools view --types snps '
         '--output-type v --output-file /dev/stdout {input.vcf} | '
         'bgzip -c /dev/stdin > {output.snps}'
         ' && '
         'bcftools view --types indels '
         '--output-type v --output-file /dev/stdout {input.vcf} | '
         'bgzip -c /dev/stdin > {output.indels}'


rule split_callset_by_genotype:
    input:
        vcf = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/20-split-vtype',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.{var_type}.vcf.bgz'
            ),
        tbi = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/20-split-vtype',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.{var_type}.vcf.bgz.tbi'
            ),
    output:
        hom = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/30-split-gtype',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.{var_type}.hom.vcf.bgz'
            ),
        het = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/30-split-gtype',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.{var_type}.het.vcf.bgz'
            ),
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
         'bcftools view --genotype hom '
         '--output-type v --output-file /dev/stdout {input.vcf} | '
         'bgzip -c /dev/stdin > {output.hom}'
         ' && '
         'bcftools view --genotype het '
         '--output-type v --output-file /dev/stdout {input.vcf} | '
         'bgzip -c /dev/stdin > {output.het}'


rule compute_vcf_stats:
    input:
        'output/evaluation/qv_estimation/variant_calls/{step}/{filename}.vcf.bgz',
        'output/evaluation/qv_estimation/variant_calls/{step}/{filename}.vcf.bgz.tbi',
    output:
        'output/evaluation/qv_estimation/variant_stats/{step}/{filename}.stats'
    conda:
         '../../environment/conda/conda_biotools.yml'
    shell:
        'bcftools stats {input[0]} > {output}'
