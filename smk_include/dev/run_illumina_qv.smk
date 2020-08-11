
QVEST_CONFIG = {
    'discard_flag': 2816,  # not primary OR qc fail OR supplementary
    'cov_threshold': '{mean} + 3 * {stddev}',  # avoid that freebayes gets tangled in ultra high cov regions
    'genome_size': int(3.1e9),  # this is roughly the sequence length of the HGSVC2 reference
    'freebayes_timeout_sec': 14400,
    'variant_min_depth': 3,
    'ref_assembly': 'GRCh38_HGSVC2_noalt',
    'skip_short_read_sources': ['PRJEB3381', 'PRJEB9396', 'HPG', 'PTG'],
    'PAV_version': 3,
    'reciprocal_overlap': [50, 25]
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


def illumina_qv_determine_targets(wildcards):
    """
    NA19239_hgsvc_pbsq2-clr_1000-flye.h2-un.arrow-p1.fasta
    """
    module_outputs = {
        'bam_stats': os.path.join('output/alignments/short_to_phased_assembly',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.mdup.stats'),
        'raw_calls': os.path.join('output/evaluation/qv_estimation/variant_stats/00-raw',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.stats'),
        'hom_snps': os.path.join('output/evaluation/qv_estimation/variant_stats/30-split-gtype',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.snvs.hom.stats'),
        'hom_ins': os.path.join('output/evaluation/qv_estimation/variant_stats/30-split-gtype',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.indels.hom.stats'),
        'het_snps': os.path.join('output/evaluation/qv_estimation/variant_stats/30-split-gtype',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.snvs.het.stats'),
        'het_ins': os.path.join('output/evaluation/qv_estimation/variant_stats/30-split-gtype',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.indels.het.stats'),
        'summary': os.path.join('output/evaluation/qv_estimation/variant_calls/70-summary-{known_ref}',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.summary.tsv'),
        'qv_est': os.path.join('output/evaluation/qv_estimation/illumina_calls_{known_ref}',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.qv.tsv'),
        'pav_isect_in': os.path.join('output/evaluation/qv_estimation/variant_calls/65-PAV-{known_ref}',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.{var_type}.{genotype}.PAV-dropped-v{pav_freeze}.r{recip_ovl}.in-hc.bed'),
        'pav_isect_out': os.path.join('output/evaluation/qv_estimation/variant_calls/65-PAV-{known_ref}',
                        '{sample}_{readset}_map-to_{assembly}.{hap}.{polisher}.{var_type}.{genotype}.PAV-dropped-v{pav_freeze}.r{recip_ovl}.out-hc.bed')
    }

    fix_wildcards = {
        'known_ref': QVEST_CONFIG['ref_assembly'],
        'pav_freeze': QVEST_CONFIG['PAV_version']
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
                if target.endswith('stats'):
                    continue
                if 'PAV' in target:
                    for rovl in QVEST_CONFIG['reciprocal_overlap']:
                        for vt, gt in [('ins', 'hom'), ('ins', 'het'), ('dels', 'hom'), ('dels', 'het')]:
                            tmp['recip_ovl'] = rovl
                            tmp['var_type'] = vt
                            tmp['genotype'] = gt
                            fmt_target = target.format(**tmp)
                            compute_results.add(fmt_target)
                else:
                    fmt_target = target.format(**tmp)
                    compute_results.add(fmt_target)
    
    return sorted(compute_results)


localrules: master_qv_estimate,
            compute_coverage_statistics,
            convert_faidx_to_regions


rule master_qv_estimate:
    input:
        illumina_qv_determine_targets


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
        mem_per_cpu_mb = lambda wildcards, attempt: int((65536 + 24576 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 65536 + 24576 * attempt,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
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
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
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
            'output/evaluation/qv_estimation/variant_calls/10-qfilter',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.q10.vcf.bgz'
            ),
    log:
        vcf = os.path.join(
            'log', 'output/evaluation/qv_estimation/variant_calls/10-qfilter',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.q10.log'
            ),
    wildcard_constraints:
        individual = '[A-Z0-9]+'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'bcftools filter --output-type u --output /dev/stdout '
        ' --include "QUAL>=10" {input.vcf} 2> {log} | '
        ' bcftools norm --fasta-ref {input.ref} --output /dev/stdout --output-type v /dev/stdin | '
        ' bgzip -c /dev/stdin > {output}'


rule split_callset_by_variant_type:
    input:
        vcf = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/10-qfilter',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.q10.vcf.bgz'
            ),
        tbi = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/10-qfilter',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.q10.vcf.bgz.tbi'
            ),
    output:
        snps = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/20-split-vtype',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.snvs.vcf.bgz'
            ),
        indels = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/20-split-vtype',
            '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.indels.vcf.bgz'
            ),
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
         'bcftools view --types snps --min-alleles 2 --max-alleles 2 '
         '--output-type v --output-file /dev/stdout {input.vcf} | '
         'bgzip -c /dev/stdin > {output.snps}'
         ' && '
         'bcftools view --types indels --min-alleles 2 --max-alleles 2 '
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
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
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
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'bcftools stats {input[0]} > {output}'


####################################################
# Below this point: lift variants to hg38 for
# evaluation restricted to high confidence regions
####################################################


rule build_contig_to_reference_paf_alignment:
    """
    Usage note: most paftools do not support input generated with
    the "--eqx" option; the below call is thus different from the
    pipeline version producing SAM/BAM output
    """
    input:
        ref = 'references/assemblies/{known_ref}.fasta',
        assm = os.path.join('output/evaluation/phased_assemblies',
                '{individual}_{assembly}.{hap}.{polisher}.fasta')
    output:
        'output/alignments/contigs_to_reference/evaluation/{individual}_{assembly}.{hap}.{polisher}_map-to_{known_ref}.paf'
    log:
        'log/output/alignments/contigs_to_reference/evaluation/{individual}_{assembly}.{hap}.{polisher}_map-to_{known_ref}.minimap.log'
    benchmark:
        'run/output/alignments/contigs_to_reference/evaluation/{individual}_{assembly}.{hap}.{polisher}_map-to_{known_ref}.minimap' + '.t{}.rsrc'.format(config['num_cpu_high'])
    conda:
        '../../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        individual = '[A-Z0-9]+'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((32768 + 8192 * attempt) / config['num_cpu_high']),
        mem_total_mb = lambda wildcards, attempt: 32768 + 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'minimap2 -t {threads} -cx asm20 --cs '
        '--secondary=no -Y -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 '
        '{input.ref} {input.assm} > {output} 2> {log}'


rule dump_vcf_snv_to_bed:
    input:
        snvs = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/30-split-gtype',
            '{callset}.snvs.{genotype}.vcf.bgz'),
    output:
        'output/evaluation/qv_estimation/variant_calls/40-dump-bed/{callset}.snvs.{genotype}.vcf.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
    shell:
        'bgzip -d -c < {input} | vcf2bed --do-not-split --snvs > {output}'


rule dump_vcf_insertions_to_bed:
    input:
        indels = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/30-split-gtype',
            '{callset}.indels.{genotype}.vcf.bgz'),
    output:
        'output/evaluation/qv_estimation/variant_calls/40-dump-bed/{callset}.ins.{genotype}.vcf.bed'
    conda:
         '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
    shell:
         'bgzip -d -c < {input} | vcf2bed --do-not-split --insertions > {output}'


rule dump_vcf_deletions_to_bed:
    input:
        indels = os.path.join(
            'output/evaluation/qv_estimation/variant_calls/30-split-gtype',
            '{callset}.indels.{genotype}.vcf.bgz'),
    output:
        'output/evaluation/qv_estimation/variant_calls/40-dump-bed/{callset}.dels.{genotype}.vcf.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
    shell:
        'bgzip -d -c < {input} | vcf2bed --do-not-split --deletions > {output}'


rule lift_call_sets_to_reference:
    input:
        bed = 'output/evaluation/qv_estimation/variant_calls/40-dump-bed/{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.{var_type}.{genotype}.vcf.bed',
        paf = 'output/alignments/contigs_to_reference/evaluation/{individual}_{assembly}.{hap}.{polisher}_map-to_{known_ref}.paf'
    output:
        'output/evaluation/qv_estimation/variant_calls/50-lifted-{known_ref}/{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.{var_type}.{genotype}.vcf.bed'
    wildcard_constraints:
        individual = '[A-Z0-9]+'
    conda:
        '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'paftools.js liftover -l 10000 {input.paf} {input.bed} > {output}'


rule restrict_calls_to_high_conf_regions:
    input:
        calls = os.path.join('output/evaluation/qv_estimation/variant_calls',
                    '50-lifted-{known_ref}',
                    '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.{var_type}.{genotype}.vcf.bed'),
        regions = 'references/annotation/hg38_giab_highconf.bed'
    output:
        hc_in = os.path.join('output/evaluation/qv_estimation/variant_calls',
                    '60-highconf-{known_ref}',
                    '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.{var_type}.{genotype}.in-hc.bed'
                ),
        hc_out = os.path.join('output/evaluation/qv_estimation/variant_calls',
                    '60-highconf-{known_ref}',
                    '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.{var_type}.{genotype}.out-hc.bed')
    conda:
         '../../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        'bedtools intersect -u -a {input.calls} -b {input.regions} > {output.hc_in}'
        ' && '
        'bedtools intersect -v -a {input.calls} -b {input.regions} > {output.hc_out}'


rule intersect_illumina_assembly_variants:
    input:
        pav = 'references/annotation/PAV_sv-insdel-dropped_v{freeze_version}.4c.bed',
        hc_in = rules.restrict_calls_to_high_conf_regions.output.hc_in,
        hc_out = rules.restrict_calls_to_high_conf_regions.output.hc_out
    output:
        isec_hc_in = os.path.join('output/evaluation/qv_estimation/variant_calls',
                        '65-PAV-{known_ref}',
                        '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.{var_type}.{genotype}.PAV-dropped-v{freeze_version}.r{recip_ovl}.in-hc.bed'),
        isec_hc_out = os.path.join('output/evaluation/qv_estimation/variant_calls',
                        '65-PAV-{known_ref}',
                        '{individual}_{library_id}_short_map-to_{assembly}.{hap}.{polisher}.{var_type}.{genotype}.PAV-dropped-v{freeze_version}.r{recip_ovl}.out-hc.bed'),
    conda:
        '../../environment/conda/conda_biotools.yml'
    params:
        recip_ovl = lambda wildcards: str(round(int(wildcards.recip_ovl) / 100, 2))
    shell:
        'bedtools intersect -wao -f {params.recip_ovl} -r -a {input.pav} -b {input.hc_in} > {output.isec_hc_in}'
        ' && '
        'bedtools intersect -wao -f {params.recip_ovl} -r -a {input.pav} -b {input.hc_out} > {output.isec_hc_out}'


def parse_full_vcf_bed_line(line):
    import re
    seq_id, start, end, _, quality, ref_allele, alt_allele, _, attributes, sample_format, sample_values = line.strip().split()
    assert sample_format.split(':')[1] == 'DP', 'Malformed VCF line in BED: {}'.format(line.strip())
    read_depth = int(sample_values.split(':')[1])
    mobj = re.search('LEN=[0-9]+;', attributes)
    if mobj is None:
        raise ValueError('No variant length in variant attributes: {}'.format(line.strip()))
    variant_length = mobj.group(0)
    variant_length = int(variant_length.strip('LEN=;'))
    return (seq_id, int(start), int(end)), (float(quality), read_depth, variant_length, ref_allele, alt_allele)


def parse_lifted_vcf_bed_line(line):

    ref_seq, ref_start, ref_end, source_location, _, _ = line.strip().split()
    seq_id, start, end = source_location.rsplit('_', 2)
    ref_start = int(ref_start)
    ref_end = int(ref_end)
    try:
        start = int(start)
        end = int(end)
    except ValueError:
        # sometimes weird suffix "_t5" or "_t3" at source location
        source_location = source_location.rsplit('_', 1)[0]
        seq_id, start, end = source_location.rsplit('_', 2)
        start = int(start)
        end = int(end)
        #raise ValueError('Cannot process line:\n\n{}\n\n'.format(line.strip()))
    return (seq_id, start, end), (ref_seq, ref_start, ref_end)


def build_complete_variant_table(variant_bed_files):

    import pandas as pd

    final_concat = []

    line_parser = parse_full_vcf_bed_line

    for bedfile in variant_bed_files:
        rows = []
        index = []
        fname = os.path.basename(bedfile)
        idx_var, idx_gt = fname.split('.')[-4:-2]
        if len(idx_var) > 3:
            idx_var = idx_var.rstrip('s')
        with open(bedfile, 'r') as bed:
            for line in bed:
                line_idx, line_values = line_parser(line)
                rows.append(line_values)
                index.append((idx_var, idx_gt, *line_idx))
        bed_index = pd.MultiIndex.from_tuples(
                        index,
                        names=['variant', 'genotype', 'seq', 'start', 'end']
                    )
        bed_df = pd.DataFrame.from_records(
                    rows,
                    index=bed_index,
                    columns=['quality', 'depth', 'variant_length', 'ref_allele', 'alt_allele']
                )
        final_concat.append(bed_df)
    
    return pd.concat(final_concat, axis=0)


def build_lifted_variant_table(lifted_bed_files):

    import pandas as pd

    final_concat = []

    line_parser = parse_lifted_vcf_bed_line

    for bedfile in lifted_bed_files:
        rows = []
        index = []
        fname = os.path.basename(bedfile)
        idx_var, idx_gt = fname.split('.')[-4:-2]
        if len(idx_var) > 3:
            idx_var = idx_var.rstrip('s')
        is_hc = 1 if fname.endswith('in-hc.bed') else 0
        with open(bedfile, 'r') as bed:
            for line in bed:
                line_idx, line_values = line_parser(line)
                rows.append((*line_values, 1, is_hc))
                index.append((idx_var, idx_gt, *line_idx))
        bed_index = pd.MultiIndex.from_tuples(
                        index,
                        names=['variant', 'genotype', 'seq', 'start', 'end']
                    )
        bed_df = pd.DataFrame.from_records(
                    rows,
                    index=bed_index,
                    columns=['ref_seq', 'ref_start', 'ref_end', 'is_lifted', 'is_highconf']
                )
        final_concat.append(bed_df)
    
    return pd.concat(final_concat, axis=0)


rule summarize_variant_calls:
    input:
        hap_assm = expand('output/evaluation/qv_estimation/variant_calls/40-dump-bed/{{callset}}.{var_type}.{genotype}.vcf.bed',
                            var_type=['snvs', 'ins', 'dels'],
                            genotype=['hom', 'het']),
        ref_assm = expand('output/evaluation/qv_estimation/variant_calls/60-highconf-{{known_ref}}/{{callset}}.{var_type}.{genotype}.{location}.bed',
                            var_type=['snvs', 'ins', 'dels'],
                            genotype=['hom', 'het'],
                            location=['in-hc', 'out-hc'])
    output:
        'output/evaluation/qv_estimation/variant_calls/70-summary-{known_ref}/{callset}.summary.tsv'
    benchmark:
        'run/output/evaluation/qv_estimation/variant_calls/70-summary-{known_ref}/{callset}.summary.rsrc'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
    run:
        import pandas as pd

        full_set = build_complete_variant_table(input.hap_assm)

        lift_set = build_lifted_variant_table(input.ref_assm)

        full_set = full_set.join(lift_set, how='outer')

        full_set.loc[full_set['ref_seq'].isna(), 'ref_seq'] = 'none'
        full_set.loc[full_set['ref_start'].isna(), 'ref_start'] = 0
        full_set.loc[full_set['ref_end'].isna(), 'ref_end'] = 0
        full_set.loc[full_set['is_lifted'].isna(), 'is_lifted'] = 0
        full_set.loc[full_set['is_highconf'].isna(), 'is_highconf'] = 0

        full_set['ref_start'] = full_set['ref_start'].astype('int64')
        full_set['ref_end'] = full_set['ref_end'].astype('int64')
        full_set['is_lifted'] = full_set['is_lifted'].astype('bool')
        full_set['is_highconf'] = full_set['is_highconf'].astype('bool')

        full_set.to_csv(output[0], sep='\t', na_rep='n/a')
    # END OF RUN BLOCK


def comp_qv(num_error_bp, bp_ref=QVEST_CONFIG['genome_size']):
    import math
    p = (num_error_bp / (bp_ref * 2))
    try:
        q = -10 * math.log10(p)
    except ValueError:
        return 100
    return int(round(q, 0))


rule compute_illumina_qv_estimate:
    input:
        'output/evaluation/qv_estimation/variant_calls/70-summary-{known_ref}/{callset}.summary.tsv'
    output:
        'output/evaluation/qv_estimation/illumina_calls_{known_ref}/{callset}.qv.tsv'
    benchmark:
        'run/output/evaluation/qv_estimation/illumina_calls_{known_ref}/{callset}.qv.rsrc'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
    run:
        import pandas as pd
        import collections as col

        index_cols = ['variant', 'genotype', 'seq', 'start', 'end']

        df = pd.read_csv(
            input[0],
            sep='\t',
            index_col=False
            )
        
        # why not setting index when reading table?
        # avoid annoying FutureWarning
        # stackoverflow.com/questions/40659212/futurewarning-elementwise-comparison-failed-returning-scalar-but-in-the-futur/46721064#46721064

        df.set_index(index_cols, inplace=True)
        
        stats = col.OrderedDict()

        stats['variants_total_count'] = df.shape[0]
        stats['variants_total_length'] = int(df['variant_length'].sum())

        df = df.loc[df['depth'] > (QVEST_CONFIG['variant_min_depth'] - 1), :].copy()

        stats['variants_total_DPgeq{}_count'.format(QVEST_CONFIG['variant_min_depth'])] = df.shape[0]
        stats['variants_total_DPgeq{}_length'.format(QVEST_CONFIG['variant_min_depth'])] = int(df['variant_length'].sum())

        for variant in ['snv', 'ins', 'del']:
            subset = df.xs(variant, level='variant')
            stats['{}_total_count'.format(variant)] = subset.shape[0]
            stats['{}_total_length'.format(variant)] = int(subset['variant_length'].sum())
            
            is_lifted = subset['is_lifted']
            
            stats['{}_lifted_count'.format(variant)] = subset.loc[is_lifted, :].shape[0]
            stats['{}_lifted_length'.format(variant)] = int(subset.loc[is_lifted, 'variant_length'].sum())
            
            stats['{}_unlifted_count'.format(variant)] = subset.loc[~is_lifted, :].shape[0]
            stats['{}_unlifted_length'.format(variant)] = int(subset.loc[~is_lifted, 'variant_length'].sum())
            
            pct = stats['{}_lifted_count'.format(variant)] / stats['{}_total_count'.format(variant)]
            stats['{}_lifted_pct'.format(variant)] = round(pct * 100, 2)
            
            is_hc = subset['is_highconf']
            
            stats['{}_highconf_count'.format(variant)] = subset.loc[is_hc, :].shape[0]
            stats['{}_highconf_length'.format(variant)] = int(subset.loc[is_hc, 'variant_length'].sum())
            pct = stats['{}_highconf_count'.format(variant)] / stats['{}_lifted_count'.format(variant)]
            stats['{}_highconf_pct'.format(variant)] = round(pct * 100, 2)
                
            for genotype in ['hom', 'het']:
                subsubset = subset.xs(genotype, level='genotype')
                stats['{}_{}_count'.format(variant, genotype)] = subsubset.shape[0]
                stats['{}_{}_length'.format(variant, genotype)] = int(subsubset['variant_length'].sum())
                
                is_lifted = subsubset['is_lifted']
            
                stats['{}_{}_lifted_count'.format(variant, genotype)] = subsubset.loc[is_lifted, :].shape[0]
                stats['{}_{}_lifted_length'.format(variant, genotype)] = int(subsubset.loc[is_lifted, 'variant_length'].sum())
                
                stats['{}_{}_unlifted_count'.format(variant, genotype)] = subsubset.loc[~is_lifted, :].shape[0]
                stats['{}_{}_unlifted_length'.format(variant, genotype)] = int(subsubset.loc[~is_lifted, 'variant_length'].sum())
                
                pct = stats['{}_{}_lifted_count'.format(variant, genotype)] / stats['{}_{}_count'.format(variant, genotype)]
                stats['{}_{}_lifted_pct'.format(variant, genotype)] = round(pct * 100, 2)

                is_hc = subsubset['is_highconf']

                stats['{}_{}_highconf_count'.format(variant, genotype)] = subsubset.loc[is_hc, :].shape[0]
                stats['{}_{}_highconf_length'.format(variant, genotype)] = int(subsubset.loc[is_hc, 'variant_length'].sum())
                pct = stats['{}_{}_highconf_count'.format(variant, genotype)] / stats['{}_{}_lifted_count'.format(variant, genotype)]
                stats['{}_{}_highconf_pct'.format(variant, genotype)] = round(pct * 100, 2)
                
            
        subset = df.xs('hom', level='genotype')
        stats['hom_highconf_count'] = subset.loc[subset['is_highconf'], :].shape[0]
        stats['hom_highconf_length'] = int(subset.loc[subset['is_highconf'], 'variant_length'].sum())

        stats['hom_unlifted_count'] = subset.loc[~subset['is_lifted'], :].shape[0]
        stats['hom_unlifted_length'] = int(subset.loc[~subset['is_lifted'], 'variant_length'].sum())

        qv_labels = [
            'QV_all_hom_highconf',
            'QV_all_hom_highconf|unlifted',
            'QV_snv_hom_highconf',
            'QV_snv_hom_highconf|unlifted',
            'QV_indels_hom_highconf',
            'QV_indels_hom_highconf|unlifted'
        ]

        qv_stat_selectors = [
            ('hom_highconf_length',),
            ('hom_highconf_length', 'hom_unlifted_length'),
            ('snv_hom_highconf_length',),
            ('snv_hom_highconf_length', 'snv_hom_unlifted_length'),
            ('ins_hom_highconf_length', 'del_hom_highconf_length'),
            ('ins_hom_highconf_length', 'ins_hom_unlifted_length', 'del_hom_highconf_length', 'del_hom_unlifted_length')
        ]

        for label, stat_select in zip(qv_labels, qv_stat_selectors):
            num_bp = sum(stats[s] for s in stat_select)
            stats[label] = comp_qv(num_bp)
            stats[label + '_bp'] = num_bp
        
        with open(output[0], 'w') as dump:
            _ = dump.write('\n'.join(['{}\t{}'.format(k, v) for k, v in stats.items()]) + '\n')

    # END OF RUN BLOCK
