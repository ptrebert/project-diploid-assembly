
localrules: master_variant_calling,
            write_intermediate_vcf_splits,
            write_final_vcf_splits


rule master_variant_calling:
    input:
        []


rule compute_position_coverage:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        # require FASTA index here because bedtools fails otherwise
        ref_idx = 'output/reference_assembly/{folder_path}/{reference}.fasta.fai',
        read_ref_aln = 'output/alignments/reads_to_reference/{folder_path}/{vc_reads}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/{folder_path}/{vc_reads}_map-to_{reference}.psort.sam.bam.bai'
    output:
        'output/alignments/reads_to_reference/{folder_path}/aux_files/{vc_reads}_map-to_{reference}.pos-cov.txt.gz'
    benchmark:
        'rsrc/output/alignments/reads_to_reference/{folder_path}/aux_files/{vc_reads}_map-to_{reference}.pos-cov.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: 2
    resources:
        mem_total_mb = lambda wildcards, attempt: 2048 + 2048 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: (2048 + 2048 * attempt) // 2
    shell:
        'bedtools genomecov -d -ibam {input.read_ref_aln} | gzip > {output}'


rule compute_uniform_coverage_regions:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        seq_info = 'output/reference_assembly/{folder_path}/{reference}/sequences/{sequence}.seq',
        pos_cov = rules.compute_position_coverage.output[0]
    output:
        'output/alignments/reads_to_reference/{folder_path}/aux_files/{vc_reads}_map-to_{reference}/{sequence}.unicov.regions'
    log:
        'log/output/alignments/reads_to_reference/{folder_path}/aux_files/{vc_reads}_map-to_{reference}/{sequence}.unicov.log'
    benchmark:
        'rsrc/output/alignments/reads_to_reference/{folder_path}/aux_files/{vc_reads}_map-to_{reference}/{sequence}.unicov.rsrc'
    conda:
        '../environment/conda/conda_pyscript.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 + 4096 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 + 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        num_regions = 128,
        script_exec = lambda wildcards: find_script_path('np_cov_to_regions.py')
    shell:
        'zgrep "{wildcards.sequence}\s" {input.pos_cov} | {params.script_exec} --debug '
            ' --seq-info {input.seq_info} --num-regions {params.num_regions} --output {output}'
            ' &> {log}'


rule call_variants_freebayes_parallel:
    """
    NB: This works on separate region files, and in parallel on each file.
        In other words, for a substantial speed up, this needs to be
        distributed across a cluster or, at least, this job should run
        on a large compute node

    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        read_ref_aln = 'output/alignments/reads_to_reference/clustered/{sseq_reads}/{vc_reads}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/clustered/{sseq_reads}/{vc_reads}_map-to_{reference}.psort.sam.bam.bai',
        reference = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta',
        ref_idx = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta.fai',
        ref_regions = 'output/alignments/reads_to_reference/clustered/{sseq_reads}/aux_files/{vc_reads}_map-to_{reference}/{sequence}.unicov.regions'
    output:
        'output/variant_calls/freebayes/{reference}/{sseq_reads}/processing/10-norm/splits/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/freebayes/{reference}/{sseq_reads}/processing/10-norm/splits/{vc_reads}.{sequence}.log'
    benchmark:
        'rsrc/output/variant_calls/freebayes/{reference}/{sseq_reads}/processing/10-norm/splits/{vc_reads}.{sequence}.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((16384 if attempt <= 1 else 16384 + attempt**attempt * 16384) / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 16384 if attempt <= 1 else 16384 + attempt**attempt * 16384,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    params:
        timeout = config['freebayes_timeout_sec'],
        script_exec = lambda wildcards: find_script_path('fb-parallel-timeout.sh')
    shell:
        '{params.script_exec} {input.ref_regions} {threads} {params.timeout} {log}'
            ' --use-best-n-alleles 4 --strict-vcf -f {input.reference} {input.read_ref_aln} > {output}'


rule call_variants_longshot:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        read_ref_aln = 'output/alignments/reads_to_reference/clustered/{sseq_reads}/{vc_reads}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/clustered/{sseq_reads}/{vc_reads}_map-to_{reference}.psort.sam.bam.bai',
        reference = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta',
        ref_idx = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta.fai',
        seq_info = 'output/reference_assembly/clustered/{sseq_reads}/{reference}/sequences/{sequence}.seq',
    output:
        'output/variant_calls/longshot/{reference}/{sseq_reads}/processing/00-raw/splits/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/longshot/{reference}/{sseq_reads}/processing/00-raw/splits/{vc_reads}.{sequence}.log'
    benchmark:
        'rsrc/output/variant_calls/longshot/{reference}/{sseq_reads}/processing/00-raw/splits/{vc_reads}.{sequence}.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    params:
        individual = lambda wildcards: wildcards.vc_reads.split('_')[0]
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 + 4096 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 + 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 6 * attempt
    shell:
        'longshot --no_haps --force_overwrite --auto_max_cov --bam {input.read_ref_aln} '
            ' --ref {input.reference} --region {wildcards.sequence}'
             ' --sample_id {params.individual} --out {output} &> {log}'


rule normalize_longshot_vcf:
    """
    NB: Currently, longshot does not adhere to the VCF specification
        and uses a "." as decimal mark for the GQ format field. WhatsHap
        cannot handle this. This step simply rounds the values to full
        integers.

    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        'output/variant_calls/longshot/{reference}/{sseq_reads}/processing/00-raw/splits/{vc_reads}.{sequence}.vcf'
    output:
        'output/variant_calls/longshot/{reference}/{sseq_reads}/processing/10-norm/splits/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/longshot/{reference}/{sseq_reads}/processing/10-norm/splits/{vc_reads}.{sequence}.log'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 + 1024 * attempt,
        mem_total_mb = lambda wildcards, attempt: 1024 + 1024 * attempt,
    run:
        import io
        with open(log[0], 'w') as logfile:
            _ = logfile.write('Processing longshot VCF from path: {}\n'.format(input[0]))

            header_lines = 0
            record_lines = 0
            filtered_dense = 0
            filtered_depth = 0

            out_buffer = io.StringIO()
            with open(input[0], 'r') as vcf_in:
                for line in vcf_in:
                    if line.startswith('##FORMAT=<ID=GQ,Number=1,Type=Float,'):
                        header_lines += 1
                        _ = logfile.write('Found GQ format field of type float\n')
                        out_buffer.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
                    elif not line.startswith('#'):
                        record_lines += 1
                        cols = line.strip().split('\t')
                        # filter out likely false positives marked as
                        # "dense cluster (dn)" in the FILTER field
                        if 'dn' in cols[6]:
                            filtered_dense += 1
                            continue
                        elif 'dp' in cols[6]:
                            # should not occur too frequently because
                            # of auto_max_cov parameter
                            filtered_depth += 1
                            continue
                        else:
                            # no other longshot qualifier known atm
                            pass
                        genotype_format = cols[-2].split(':')
                        if not genotype_format[1] == 'GQ':
                            _ = logfile.write('ERROR - approximate line number: {}\n'.format(header_lines + record_lines))
                            _ = logfile.write('Unexpected FORMAT composition (GQ not at position 1): {}\n'.format(line.strip()))
                            _ = logfile.write('ERROR - abort\n')
                            raise ValueError('Unexpected FORMAT field, check log file: {}\n'.format(log[0]))
                        genotype_values = cols[-1].split(':')
                        genotype_values[1] = str(int(round(float(genotype_values[1]), 0)))
                        genotype_values = ':'.join(genotype_values)
                        new_line = '\t'.join(cols[:-1]) + '\t' + genotype_values
                        out_buffer.write(new_line + '\n')
                    else:
                        header_lines += 1
                        out_buffer.write(line)
            _ = logfile.write('All lines of input VCF buffered (header: {} / records: {})\n'.format(header_lines, record_lines))
            _ = logfile.write('Header lines: {}\n'.format(header_lines))
            _ = logfile.write('Records total: {}\n'.format(record_lines))
            _ = logfile.write('SNVs in dense clusters removed: {}\n'.format(filtered_dense))
            _ = logfile.write('SNVs exceeding max depth removed: {}\n'.format(filtered_depth))
            _ = logfile.write('SNV records remaining: {}\n'.format(record_lines - filtered_dense - filtered_depth))

            with open(output[0], 'w') as vcf_out:
                _ = vcf_out.write(out_buffer.getvalue())

            _ = logfile.write('Normalized longshot VCF record dumped to file: {}\n'.format(output[0]))
            _ = logfile.write('Done')
    # end of run block


rule call_variants_deepvariant:
    """
    Dangerous to run several DeepVariant jobs on same server, see here:
    github.com/google/deepvariant/issues/242
    Suggested workaround:
    --intermediate_results_dir="/tmp/deepvariant_tmp_output/chr1"
    So, keep this as part of the call to schedule more DeepVariant jobs in parallel

    For DeepVariant >= 1.0:
    use of phase information for variant calling requires the HP tag to be set in the
    BAM file. This is generally not the case for the assembly pipeline (= phasing happens
    after variant calling).
    """
    input:
        container = 'output/container/docker/google/deepvariant_{}.sif'.format(config['deepvariant_version']),
        reference = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta',
        ref_idx = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta.fai',
        seq_info = 'output/reference_assembly/clustered/{sseq_reads}/{reference}/sequences/{sequence}.seq',
        read_ref_aln = 'output/alignments/reads_to_reference/clustered/{sseq_reads}/{vc_reads}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/clustered/{sseq_reads}/{vc_reads}_map-to_{reference}.psort.sam.bam.bai'
    output:
        vcf = 'output/variant_calls/deepvar/{reference}/{sseq_reads}/processing/10-norm/splits/{vc_reads}.{sequence}.vcf',
        gvcf = 'output/variant_calls/deepvar/{reference}/{sseq_reads}/processing/10-norm/splits/{vc_reads}.{sequence}.gvcf',
    log:
        'log/output/variant_calls/deepvar/{reference}/{sseq_reads}/processing/10-norm/splits/{vc_reads}.{sequence}.log'
    benchmark:
        os.path.join('rsrc/output/variant_calls/deepvar',
                     '{reference}/{sseq_reads}/processing/10-norm/splits',
                     '{vc_reads}.{sequence}' + '.t{}.rsrc'.format(config['num_cpu_medium']))
    envmodules:
        config['env_module_singularity']
    threads: config['num_cpu_medium']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: int((16384 + (max(attempt - 1, 0) * 16384)) / config['num_cpu_medium']),
        mem_total_mb = lambda wildcards, attempt: 16384 + (max(attempt - 1, 0) * 16384),
        runtime_hrs = lambda wildcards, attempt: attempt * attempt * attempt
    params:
        bind_folder = lambda wildcards: os.getcwd(),
        temp_dir = lambda wildcards: os.path.join('/tmp', 'deepvariant', wildcards.reference, wildcards.sseq_reads, wildcards.vc_reads, wildcards.sequence),
        use_hap_info = '--nouse_hp_information' if config['git_commit_version'] > 12 else '',
        singularity = '' if not config.get('env_module_singularity', False) else 'module load {} ; '.format(config['env_module_singularity'])
    shell:
        '{params.singularity}'
        'singularity run --bind {params.bind_folder}:/wd {input.container} /opt/deepvariant/bin/run_deepvariant '
            ' --model_type=PACBIO  --ref=/wd/{input.reference} --reads=/wd/{input.read_ref_aln} '
            ' --regions "{wildcards.sequence}" --output_vcf=/wd/{output.vcf} --output_gvcf=/wd/{output.gvcf} '
            ' {params.use_hap_info} --novcf_stats_report '
            ' --intermediate_results_dir="{params.temp_dir}" --num_shards={threads} &> {log} ; '
        'rm -rfd {params.temp_dir} ; '


rule filter_variant_calls_quality_biallelic_snps:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/10-norm/splits/{vc_reads}.{sequence}.vcf.bgz',
        tbi = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/10-norm/splits/{vc_reads}.{sequence}.vcf.bgz.tbi',
        reference = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta'
    output:
        'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/splits/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/splits/{vc_reads}.{sequence}.log'
    conda:
        '../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 + 1024 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 + 1024 * attempt
    shell:
        'bcftools filter --include "QUAL>={wildcards.qual}" {input.vcf} | '
            'bcftools view -c 1 --types snps -m 2 -M 2 | '
            'bcftools norm --fasta-ref {input.reference} --output {output} --output-type v &> {log}'


rule whatshap_regenotype_variant_calls:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/splits/{vc_reads}.{sequence}.vcf',
        read_ref_aln = 'output/alignments/reads_to_reference/clustered/{sseq_reads}/{vc_reads}_map-to_{reference}.psort.sam.bam',
        aln_idx = 'output/alignments/reads_to_reference/clustered/{sseq_reads}/{vc_reads}_map-to_{reference}.psort.sam.bam.bai',
        reference = 'output/reference_assembly/clustered/{sseq_reads}/{reference}.fasta'
    output:
       'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/30-regenotype/splits/{vc_reads}.{sequence}.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/30-regenotype/splits/{vc_reads}.{sequence}.log'
    benchmark:
        'rsrc/output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/30-regenotype/splits/{vc_reads}.{sequence}.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 2048 + 1024 * attempt,
        mem_total_mb = lambda wildcards, attempt: 2048 + 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: 3 if attempt <= 1 else 12 * attempt
    shell:
        'whatshap genotype --chromosome {wildcards.sequence} --reference {input.reference} --output {output} {input.vcf} {input.read_ref_aln} &> {log}'


rule filter_retyped_by_genotype_quality:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/30-regenotype/splits/{vc_reads}.{sequence}.vcf'
    output:
        temp('output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/35-filter-GQ{gq}/splits/{vc_reads}.{sequence}.vcf')
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/35-filter-GQ{gq}/splits/{vc_reads}.{sequence}.log'
    conda:
        '../environment/conda/conda_biotools.yml'
    shell:
        'bcftools filter --include "GQ>={wildcards.gq}" --output-type v --output {output} {input.vcf} &> {log}'


rule extract_heterozygous_variants:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        vcf_original = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/splits/{vc_reads}.{sequence}.vcf',
        vcf_retyped = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/35-filter-GQ{gq}/splits/{vc_reads}.{sequence}.vcf'
    output:
        vcf_original = temp('output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/40-extract-het-GQ{gq}/splits/{vc_reads}.{sequence}.het-only-original.vcf'),
        vcf_retyped = temp('output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/40-extract-het-GQ{gq}/splits/{vc_reads}.{sequence}.het-only-retyped.vcf'),
    conda:
        '../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 + 1024 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 + 1024 * attempt
    shell:
        'bcftools view --genotype het --output-type v --output-file {output.vcf_original} {input.vcf_original} '
        ' && '
        'bcftools view --genotype het --output-type v --output-file {output.vcf_retyped} {input.vcf_retyped} '



rule intersect_original_retyped_variant_calls:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    """
    input:
        original_vcf = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/40-extract-het-GQ{gq}/splits/{vc_reads}.{sequence}.het-only-original.vcf.bgz',
        original_tbi = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/40-extract-het-GQ{gq}/splits/{vc_reads}.{sequence}.het-only-original.vcf.bgz.tbi',
        retyped_vcf = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/40-extract-het-GQ{gq}/splits/{vc_reads}.{sequence}.het-only-retyped.vcf.bgz',
        retyped_tbi = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/40-extract-het-GQ{gq}/splits/{vc_reads}.{sequence}.het-only-retyped.vcf.bgz.tbi',
    output:
        uniq_original = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/50-intersect-GQ{gq}/splits/{vc_reads}.{sequence}/0000.vcf',
        uniq_retyped = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/50-intersect-GQ{gq}/splits/{vc_reads}.{sequence}/0001.vcf',
        shared_original = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/50-intersect-GQ{gq}/splits/{vc_reads}.{sequence}/0002.vcf',
        shared_retyped = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/50-intersect-GQ{gq}/splits/{vc_reads}.{sequence}/0003.vcf',
        desc = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/50-intersect-GQ{gq}/splits/{vc_reads}.{sequence}/README.txt',
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/50-intersect-GQ{gq}/{vc_reads}.{sequence}.isect.log'
    conda:
        '../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 + 1024 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 + 1024 * attempt
    params:
        outdir = lambda wildcards, output: os.path.dirname(output.desc)
    shell:
        'bcftools isec -p {params.outdir} {input.original_vcf} {input.retyped_vcf} &> {log}'


def collect_final_vcf_splits(wildcards, glob_collect=True, caller='snakemake'):
    """
    """
    source_path = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/50-intersect-GQ{gq}/splits/{vc_reads}.{sequence}/0002.vcf'
    glob_path = source_path.replace('.{sequence}/', '.*/')
    cluster_key = 'sequence'
    if glob_collect:
        vcf_files = check_cluster_file_completeness(wildcards, source_path, glob_path, wildcards.sseq_reads, cluster_key)
        if not vcf_files:
            raise RuntimeError('collect_final_vcf_splits: no files collected with pattern {}'.format(pattern))

    else:
        raise RuntimeError('Illegal function call: Snakemake checkpoints must not be used!')

        reference_folder = os.path.join('output/reference_assembly/clustered', wildcards.sseq_reads)
        seq_output_dir = checkpoints.create_assembly_sequence_files.get(folder_path=reference_folder,
                                                                        reference=wildcards.reference).output[0]

        checkpoint_wildcards = glob_wildcards(
            os.path.join(seq_output_dir, '{sequence}.seq')
            )

        vcf_files = expand(
            source_path,
            var_caller=wildcards.var_caller,
            reference=wildcards.reference,
            sseq_reads=wildcards.sseq_reads,
            qual=wildcards.qual,
            gq=wildcards.gq,
            vc_reads=wildcards.vc_reads,
            sequence=checkpoint_wildcards.sequence
            )

    return vcf_files


# DEPRECATED
# rule write_final_vcf_splits:
#     input:
#         vcf_splits = collect_final_vcf_splits
#     output:
#         fofn = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.fofn'
#     run:
#         import os

#         try:
#             validate_checkpoint_output(input.vcf_splits)
#             vcf_splits = input.vcf_splits
#         except (RuntimeError, ValueError) as error:
#             import sys
#             sys.stderr.write('\n{}\n'.format(str(error)))
#             vcf_splits = collect_final_vcf_splits(wildcards, glob_collect=True)


#         with open(output.fofn, 'w') as dump:
#             for file_path in sorted(vcf_splits):
#                 if not os.path.isfile(file_path):
#                     import sys
#                     sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
#                 _ = dump.write(file_path + '\n')


rule write_final_vcf_splits:
    input:
        vcf_splits = collect_final_vcf_splits
    output:
        fofn = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.fofn'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.write-fofn.log'
    run:
        import os
        import sys

        with open(log[0], 'w') as logfile:
            _ = logfile.write('write_final_vcf_splits: {}\n'.format(str(wildcards)))
            num_vcf_splits = len(input.vcf_splits)
            assert num_vcf_splits > 0, 'write_final_vcf_splits: no VCF split files as input'
            sample_name = wildcards.sseq_reads.split('_')[0]
            num_clusters = estimate_number_of_saarclusters(sample_name, wildcards.sseq_reads)
            _ = logfile.write('Number of input VCF splits: {}\n'.format(num_vcf_splits))
            _ = logfile.write('Number of expected VCF splits (clusters): {}\n'.format(num_clusters))
            _ = logfile.write('Received following VCF split input:\n{}\n'.format('\n'.join(sorted(input.vcf_splits))))
            if num_clusters != num_vcf_splits:
                _ = logfile.write('Potential error situation: unexpected number of VCF split files as input\n')
                _ = logfile.write('Checking for available VCF split files...\n')
                vcf_splits = collect_final_vcf_splits(wildcards, caller='runblock')
                if len(vcf_splits) != num_clusters:
                    _ = logfile.write('Cannot handle error: unexpected number ({} vs {}) of VCF split files\n'.format(len(vcf_splits), num_clusters))
                    _ = logfile.write('===\n{}\n'.format('\n'.join(vcf_splits)))
                    raise RuntimeError('')
            else:
                vcf_splits = input.vcf_splits

        with open(output.fofn, 'w') as fofn:
            _ = fofn.write('\n'.join(sorted(vcf_splits)))
    # END OF RUN BLOCK


rule merge_final_vcf_splits:
    input:
        fofn = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.fofn'
    output:
        'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.concat.log'
    conda:
        '../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 + 1024 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 + 1024 * attempt
    shell:
        'bcftools concat -f {input.fofn} --output {output} --output-type v &> {log}'


rule compute_final_vcf_stats:
    input:
        im_stats = 'output/statistics/variant_calls/{var_caller}/{reference}/{sseq_reads}/{vc_reads}.snv.QUAL{qual}.vcf.stats',
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.vcf.bgz',
        idx = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}_GQ{gq}/{vc_reads}.snv.vcf.bgz.tbi',
    output:
        stats = 'output/statistics/variant_calls/{var_caller}/{reference}/{sseq_reads}/{vc_reads}.snv.QUAL{qual}.GQ{gq}.vcf.stats'
    conda:
        '../environment/conda/conda_biotools.yml'
    priority: 200
    shell:
        'bcftools stats {input.vcf} > {output.stats}'


# Same again for convenience on intermediate set of SNPs

def collect_intermediate_vcf_splits(wildcards, glob_collect=True, caller='snakemake'):
    """
    """
    source_path = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/processing/20-snps-QUAL{qual}/splits/{vc_reads}.{sequence}.vcf'
    glob_path = source_path.replace('.{sequence}.', '.*.')
    cluster_key = 'sequence'
    if glob_collect:
        vcf_files = check_cluster_file_completeness(wildcards, source_path, glob_path, wildcards.sseq_reads, cluster_key)
        if not vcf_files:
            raise RuntimeError('collect_intermediate_vcf_splits: no files collected with pattern {}'.format(pattern))

    else:
        raise RuntimeError('Illegal function call: Snakemake checkpoints must not be used!')

        reference_folder = os.path.join('output/reference_assembly/clustered', wildcards.sseq_reads)
        seq_output_dir = checkpoints.create_assembly_sequence_files.get(folder_path=reference_folder,
                                                                        reference=wildcards.reference).output[0]

        checkpoint_wildcards = glob_wildcards(
            os.path.join(seq_output_dir, '{sequence}.seq')
            )

        vcf_files = expand(
            source_path,
            var_caller=wildcards.var_caller,
            reference=wildcards.reference,
            sseq_reads=wildcards.sseq_reads,
            qual=wildcards.qual,
            vc_reads=wildcards.vc_reads,
            sequence=checkpoint_wildcards.sequence
            )

    return vcf_files


# DEPRECATED
# rule write_intermediate_vcf_splits:
#     input:
#         vcf_splits = collect_intermediate_vcf_splits
#     output:
#         fofn = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}/{vc_reads}.snv.fofn'
#     run:
#         import os

#         try:
#             validate_checkpoint_output(input.vcf_splits)
#             vcf_splits = input.vcf_splits
#         except (RuntimeError, ValueError) as error:
#             import sys
#             sys.stderr.write('\n{}\n'.format(str(error)))
#             vcf_splits = collect_intermediate_vcf_splits(wildcards, glob_collect=True)

#         with open(output.fofn, 'w') as dump:
#             for file_path in sorted(vcf_splits):
#                 if not os.path.isfile(file_path):
#                     import sys
#                     sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
#                 _ = dump.write(file_path + '\n')


rule write_intermediate_vcf_splits:
    """
    2022-01-11
    Depending on when this rule is executed, it may throw an error because
    of missing VCF cluster files as a result of some clusters being dropped
    by breakpointR (see comments in integrative_phasing.smk module).
    In that case, estimate_number_of_saarclusters() would already return
    the number of clusters excluding the dropped ones. Change the code
    below to proceed in case of a number mismatch where the estimate
    is smaller than the number of input files here
    """
    input:
        vcf_splits = collect_intermediate_vcf_splits
    output:
        fofn = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}/{vc_reads}.snv.fofn'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}/{vc_reads}.write-fofn.log'
    run:
        import os
        import sys

        with open(log[0], 'w') as logfile:
            _ = logfile.write('write_intermediate_vcf_splits: {}\n'.format(str(wildcards)))
            num_vcf_splits = len(input.vcf_splits)
            assert num_vcf_splits > 0, 'write_intermediate_vcf_splits: no VCF split files as input'
            sample_name = wildcards.sseq_reads.split('_')[0]
            num_clusters = estimate_number_of_saarclusters(sample_name, wildcards.sseq_reads)
            
            _ = logfile.write('Number of input VCF splits: {}\n'.format(num_vcf_splits))
            _ = logfile.write('Number of expected VCF splits (clusters): {}\n'.format(num_clusters))
            _ = logfile.write('Received following VCF split input:\n{}\n'.format('\n'.join(sorted(input.vcf_splits))))
            if num_clusters < num_vcf_splits and ((num_clusters * 2) > num_vcf_splits):
                # if there are too many input vcf splits, this should fail nevertheless
                _ = logfile.write('Expected number of clusters is smaller, assuming clusters were dropped downstream by breakpointR\n')
                vcf_splits = input.vcf_splits
            elif num_clusters != num_vcf_splits:
                _ = logfile.write('Potential error situation: unexpected number of VCF split files as input\n')
                _ = logfile.write('Checking for available VCF split files...\n')
                vcf_splits = collect_intermediate_vcf_splits(wildcards, caller='runblock')
                if len(vcf_splits) != num_clusters:
                    _ = logfile.write('Cannot handle error: unexpected number ({} vs {}) of VCF split files\n'.format(len(vcf_splits), num_clusters))
                    _ = logfile.write('===\n{}\n'.format('\n'.join(vcf_splits)))
                    raise RuntimeError('')
            else:
                vcf_splits = input.vcf_splits

        with open(output.fofn, 'w') as fofn:
            _ = fofn.write('\n'.join(sorted(vcf_splits)))
    # END OF RUN BLOCK


rule merge_intermediate_vcf_splits:
    input:
        fofn = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}/{vc_reads}.snv.fofn'
    output:
        'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}/{vc_reads}.snv.vcf'
    log:
        'log/output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}/{vc_reads}.concat.log'
    conda:
        '../environment/conda/conda_biotools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 + 1024 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 + 1024 * attempt
    shell:
        'bcftools concat -f {input.fofn} --output {output} --output-type v &> {log}'


rule compute_intermediate_vcf_stats:
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}/{vc_reads}.snv.vcf.bgz',
        idx = 'output/variant_calls/{var_caller}/{reference}/{sseq_reads}/QUAL{qual}/{vc_reads}.snv.vcf.bgz.tbi',
    output:
        stats = 'output/statistics/variant_calls/{var_caller}/{reference}/{sseq_reads}/{vc_reads}.snv.QUAL{qual}.vcf.stats'
    conda:
        '../environment/conda/conda_biotools.yml'
    priority: 200
    shell:
        'bcftools stats {input.vcf} > {output.stats}'
