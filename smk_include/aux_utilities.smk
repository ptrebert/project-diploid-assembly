
localrules: make_chromosome_size_file, mock_index_reads


rule make_chromosome_size_file:
    input:
        'references/assemblies/{ref_genome}.fasta.fai'
    output:
        'references/assemblies/{ref_genome}.sizes'
    shell:
        "cut -f 1,2 {input} > {output}"


rule mock_index_reads:
    """
    Rule exists to create mock input for unimap
    index creation
    """
    input:
        'references/assemblies/{ref_genome}.fasta'
    output:
        temp('references/assemblies/mock_idx/{ref_genome}.index_read.fasta')
    run:
        with open(output[0], 'w') as fasta:
            _ = fasta.write('>index_read\n')
            _ = fasta.write('ACGTACGT\n')


rule create_unimap_index:
    """
    NB: index compatibility (k-mer size default: 21)
    """
    input:
        ref = 'references/assemblies/{ref_genome}.fasta',
        reads = 'references/assemblies/mock_idx/{ref_genome}.index_read.fasta'
    output:
        umi = 'references/assemblies/{ref_genome}.umi',
    log:
        'log/references/assemblies/{ref_genome}.umi.log',
    benchmark:
        'rsrc/references/assemblies/{ref_genome}.umi.rsrc',
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: 2
    resources:
        runtime_hrs = lambda wildcards, attempt: max(0, attempt - 1),
        mem_total_mb = lambda wildcards, attempt: 16384 + 16384 * attempt
    shell:
        'unimap -d {output} -x asm20 -t {threads} -o /dev/null {input.ref} {input.reads} &> {log}'
    

rule compute_md5_checksum:
    input:
        '{filepath}'
    output:
        '{filepath}.md5'
    conda:
         '../environment/conda/conda_shelltools.yml'
    shell:
        "md5sum {input} > {output}"


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
        'rsrc/{filepath}.idx-bai.t1.rsrc'
    resources:
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 16 * attempt
    conda:
         '../environment/conda/conda_biotools.yml'
    shell:
        "samtools index {input.bam}"


rule pb_bamtools_index_bam_alignment:
    input:
        pbn_bam = '{filepath}.pbn.bam'
    output:
        pbi = '{filepath}.pbn.bam.pbi'
    log:
        'log/{filepath}.create-pbi.log'
    benchmark:
        'rsrc/{filepath}.create-pbi.rsrc'
    resources:
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 16 * attempt
    conda:
         '../environment/conda/conda_pbtools.yml'
    shell:
        'pbindex {input} &> {log}'


rule pb_bam2x_dump_fastq:
    input:
        pbn_bam = 'input/bam/{pbn_sample}_{sampling}.pbn.bam',
        pbn_idx = 'input/bam/{pbn_sample}_{sampling}.pbn.bam.pbi'
    output:
        'input/fastq/{pbn_sample}_{sampling}.fastq.gz'
    log:
        'log/input/bam/{pbn_sample}_{sampling}.dump.log'
    benchmark:
        'rsrc/input/bam/{pbn_sample}_{sampling}.dump.rsrc'
    wildcard_constraints:
        pbn_sample = CONSTRAINT_ALL_PBN_INPUT_SAMPLES,
        sampling = '[0-9x]+'
    resources:
        runtime_hrs = lambda wildcards, attempt: 48 * attempt
    conda:
         '../environment/conda/conda_pbtools.yml'
    params:
        out_prefix = lambda wildcards, output: output[0].rsplit('.', 2)[0]
    shell:
        'bam2fastq -c 5 -o {params.out_prefix} {input.pbn_bam} &> {log}'


rule pb_bam2x_dump_haploid_fastq:
    input:
        pbn_bam = 'output/diploid_assembly/{folder_path}/draft/haploid_bam/{pbn_hap_reads}.{hap}.{sequence}.pbn.bam',
        pbn_idx = 'output/diploid_assembly/{folder_path}/draft/haploid_bam/{pbn_hap_reads}.{hap}.{sequence}.pbn.bam.pbi',
    output:
        'output/diploid_assembly/{folder_path}/draft/haploid_fastq/{pbn_hap_reads}.{hap}.{sequence}.fastq.gz',
    log:
        'log/output/diploid_assembly/{folder_path}/draft/haploid_fastq/{pbn_hap_reads}.{hap}.{sequence}.dump.log',
    benchmark:
        'rsrc/output/diploid_assembly/{folder_path}/draft/haploid_fastq/{pbn_hap_reads}.{hap}.{sequence}.dump.rsrc',
    wildcard_constraints:
        pbn_hap_reads = CONSTRAINT_ALL_PBN_INPUT_SAMPLES + '_[0-9]+',
    resources:
        runtime_hrs = lambda wildcards, attempt: 1 if attempt <= 1 else 2 * attempt
    conda:
         '../environment/conda/conda_pbtools.yml'
    params:
        out_prefix = lambda wildcards, output: output[0].rsplit('.', 2)[0]
    shell:
        'bam2fastq -c 5 -o {params.out_prefix} {input.pbn_bam} &> {log}'


rule samtools_index_fasta:
    input:
        '{filepath}.{fasta_ext}'
    output:
        '{filepath}.{fasta_ext}.fai'
    conda:
         '../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        fasta_ext = '(fasta|fa|fna)'
    shell:
        'samtools faidx {input}'


rule bgzip_file_copy:
    input:
        '{filepath}'
    output:
        '{filepath}.bgz'
    conda:
         '../environment/conda/conda_biotools.yml'
    shell:
        'bgzip -c {input} > {output}'


rule bcftools_index_bgzipped_file:
    input:
        '{filepath}.bgz'
    output:
        '{filepath}.bgz.tbi'
    conda:
         '../environment/conda/conda_biotools.yml'
    shell:
        'bcftools index --tbi {input}'


rule compute_genome_id:
    input:
        '{filepath}.fasta'
    output:
        '{filepath}.gid.tsv'
    conda:
         '../environment/conda/conda_pyscript.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script_exec = lambda wildcards: find_script_path('genome_id.py', 'utilities')
    shell:
        '{params.script_exec} --describe --genome {input} --output {output}'


rule create_assembly_sequence_files:
    """
    Converted from checkpoint to regular rule to get rid of
    checkpoint-related problems.
    Risk: FASTA and FASTA index may have a disconnect
    (e.g., after a --touch run) for clustered assemblies
    b/c SaaRclust-ering does not produce a fix number of clusters.
    Go extra mile here in case the requested sequence
    entry is not part of the FASTA index (make debugging easier)
    """
    input:
        fasta = '{folder_path}/{reference}.fasta',
        fai = '{folder_path}/{reference}.fasta.fai'
    output:
        '{folder_path}/{reference}/sequences/{sequence}.seq'
    run:
        import os
        import sys

        index_sequences = set()
        dump_sequence_entry = None
        with open(input.fai, 'r') as fasta_index:
            for line in fasta_index:
                sequence = line.split()[0]
                index_sequences.add(sequence)
                if sequence == wildcards.sequence:
                    dump_sequence_entry = line

        if dump_sequence_entry is None:
            sequence_records = '\n'.join(sorted(index_sequences))
            sys.stderr.write(
                'Requested sequence entry {} is not part of index file: {}\n'.format(wildcards.sequence, input.fai)
            )
            sys.stderr.write('FASTA index contains the following records:\n{}\n'.format(sequence_records))

            # now check if the FASTA contains the record
            fasta_sequences = set()
            with open(input.fasta, 'r') as fasta:
                for line in fasta:
                    if not line.startswith('>'):
                        continue
                    fasta_sequences.add(line.strip().strip('>'))
            if wildcards.sequence in fasta_sequences:
                sys.stderr.write('FASTA file contains requested sequence entry.\n')
                sys.stderr.write('Need to stop the pipeline and manually remove FASTA index file: {}\n'.format(input.fai))
            else:
                sys.stderr.write('FASTA file does not contain requested sequence entry.\n')
                sys.stderr.write('Potential fix: stop the pipeline and manually remove FASTA and FASTA index file: {} / {}\n'.format(input.fasta, input.fai))
            sequence_records = '\n'.join(sorted(fasta_sequences))
            sys.stderr.write('Found the following sequence records in FASTA file: {}\n'.format(input.fasta))
            sys.stderr.write('===\n{}\n'.format(sequence_records))
            raise ValueError('Missing sequence record: {} not in {}'.format(wildcards.sequence, input.fai))
        
        with open(output[0], 'w') as dump:
            _ = dump.write(dump_sequence_entry)
    # END OF RUN BLOCK


rule generate_bwa_index:
    input:
        reference = '{folder_path}/{file_name}.fasta'
    output:
        '{folder_path}/{file_name}/bwa_index/{file_name}.amb',
        '{folder_path}/{file_name}/bwa_index/{file_name}.ann',
        '{folder_path}/{file_name}/bwa_index/{file_name}.bwt',
        '{folder_path}/{file_name}/bwa_index/{file_name}.pac',
        '{folder_path}/{file_name}/bwa_index/{file_name}.sa'
    log:
        'log/{folder_path}/bwa_index/{file_name}.log'
    benchmark:
        'rsrc/{folder_path}/bwa_index/{file_name}.rsrc'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    conda:
         '../environment/conda/conda_biotools.yml'
    params:
        prefix = lambda wildcards, output: output[0].rsplit('.', 1)[0]
    shell:
        'bwa index -p {params.prefix} {input.reference} &> {log}'


def get_revcomp_translation_table():
    """
    Non-canonical chars (nucleotides) are returned unchanged
    when calling translate() on the sequence
    """
    revcomp_mapping = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    revcomp_mapping.update({k.lower(): v.lower() for k, v in revcomp_mapping.items()})

    revcomp_table = str.maketrans(revcomp_mapping)
    return revcomp_table


def estimate_number_of_saarclusters(sample_name, sseq_reads, return_names=False, readset=None):
    """
    Function introduced to drop all checkpoints w/ subsequent dynamic aggregation of output files.
    If number of clusters happens to match what is being produced by SaaRclust, should allow
    pipeline to run from start to end w/o interruption.

    Update 2021-12-16
    Non-consecutive naming of clusters can cause trouble. Introduce new parameter
    "return_names" with default False to not break function calls.
    """
    import os
    import sys
    import glob

    if DEBUG:
        sys.stderr.write(f'Estimating number of SaaRclusters for: {sample_name} / {sseq_reads} / {readset}\n')

    num_clusters = 0
    # 2021-12-16
    # SaaRclust (as of #4700afe) has the annoying property of not renaming clusters after
    # collapsing (?) clusters to meet the desired number of clusters;
    # as a consequence, cluster IDs are not necessarily enumerated consecutively,
    # which makes it next to impossible to detect missing information pertaining
    # to a specific cluster
    cluster_names = set()

    # added 2022-03-22
    # Triggered for downsampling project where combination of Strand-seq and LR readset
    # is not unique enough for the below glob to find the correct file (all samples have the same name).
    # Add parameter readset to disambiguate.
    if readset is None:
        # TODO this can probably be dropped
        formatter = {'sseq_reads': sseq_reads, 'sample': sample_name, 'version': config['git_commit_version']}
        cluster_fofn_path = 'output/reference_assembly/clustered/{sseq_reads}/{sample}*_scV{version}-*.cluster-ids.txt'.format(**formatter)
        # added 2021-12-24: in case breakpointR cannot process some clusters, the estimate would always be wrong, i.e.
        # need to subtract the IDs of dropped clusters
        dropped_fofn_path = 'output/reference_assembly/clustered/{sseq_reads}/{sample}*_scV{version}-*.dropped-clusters.err'.format(**formatter)
    else:
        formatter = {'sseq_reads': sseq_reads, 'readset': readset, 'version': config['git_commit_version']}
        cluster_fofn_path = 'output/reference_assembly/clustered/{sseq_reads}/{readset}_scV{version}-*.cluster-ids.txt'.format(**formatter)
        # added 2021-12-24: in case breakpointR cannot process some clusters, the estimate would always be wrong, i.e.
        # need to subtract the IDs of dropped clusters
        dropped_fofn_path = 'output/reference_assembly/clustered/{sseq_reads}/{readset}_scV{version}-*.dropped-clusters.err'.format(**formatter)

    if DEBUG:
        sys.stderr.write(f'Checking cluster fofn path: {cluster_fofn_path}\n')
    cluster_fofn_matches = glob.glob(cluster_fofn_path)
    if len(cluster_fofn_matches) == 1:
        cluster_fofn = cluster_fofn_matches[0]
        if DEBUG:
            sys.stderr.write('... glob succeeded\n')
    else:
        if len(cluster_fofn_matches) > 1:
            if WARN:
                sys.stderr.write('WARNING: ambiguous match for sequence cluster fofn file: {}\n'.format(', '.join(cluster_fofn_matches)))
        cluster_fofn = None

    if cluster_fofn is not None and os.path.isfile(cluster_fofn):
        if DEBUG:
            sys.stderr.write('Loading number of clusters from fofn: {}\n'.format(cluster_fofn))
        num_clusters_slack = int(config.get('num_cluster_slack', 1))
        default_clusters = int(config.get('num_default_clusters', 24))
        subtract_clusters = set()
        # check if dropped_clusters exists
        dropped_fofn_matches = glob.glob(dropped_fofn_path)
        if len(dropped_fofn_matches) == 1:
            dropped_fofn = dropped_fofn_matches[0]
            with open(dropped_fofn, 'r') as listing:
                subtract_clusters = set(l.strip() for l in listing if l.strip())
        # probably best case: cluster fofn has already been created
        with open(cluster_fofn, 'r') as fofn:
            for line in fofn:
                if not line.strip():
                    continue
                if 'cluster99' in line:
                    raise ValueError(f'ERROR: cluster99 in cluster FOFN file {cluster_fofn} detected. Delete this file and restart the pipeline')
                # Update 2021-12-16: in PGAS versions v14+ (later than #dedaeba), the fofn file
                # contains only the cluster names, and not the complete file paths because
                # the SaaRclust output is post-processed.
                cluster_names.add(line.strip())
        cluster_names = sorted(cluster_names - subtract_clusters)
        fofn_clusters = len(cluster_names)
        if fofn_clusters < default_clusters - num_clusters_slack:
            raise ValueError(f'ERROR: number of clusters loaded from FOFN file {cluster_fofn} is too small. Delete this file and restart the pipeline')
        num_clusters = fofn_clusters
        if DEBUG:
            sys.stderr.write('Loaded number of cluster from fofn: {}\n'.format(num_clusters))
    elif sample_name in CONSTRAINT_MALE_SAMPLE_NAMES:
        num_clusters = config.get('desired_clusters_male', config.get('desired_clusters', 0))
        if DEBUG:
            sys.stderr.write('Male sample - cluster estimate (config): {}\n'.format(num_clusters))
    elif sample_name in CONSTRAINT_FEMALE_SAMPLE_NAMES:
        num_clusters = config.get('desired_clusters_female', config.get('desired_clusters', 0))
        if DEBUG:
            sys.stderr.write('Female sample - cluster estimate (config): {}\n'.format(num_clusters))
    else:
        num_clusters = config.get('desired_clusters', 0)
        if DEBUG:
            sys.stderr.write('Generic sample - cluster estimate (config): {}\n'.format(num_clusters))
    if num_clusters == 0:
        # reasonable best guess
        num_clusters = int(config.get('num_default_clusters', 24))
        if DEBUG:
            sys.stderr.write('Cluster estimate is zero - returning best guess: {}\n'.format(num_clusters))

    if return_names:
        if not cluster_names:
            # have to work with best guess
            cluster_names = [f'cluster{i}' for i in range(1, num_clusters + 1)]
        retval = num_clusters, cluster_names
    else:
        retval = num_clusters

    return retval


def check_cluster_file_completeness(wildcards, source_path, glob_path, sseq_reads, cluster_key):
    """
    Update 2021-12-16
    Rely on cluster names returned by estimate_number_of_saarclusters
    as cluster names may not be consecutive (as-of SaaRclust #4700afe)

    # Update 2022-02-03
    There is a unresolvable mismatch between the number of clusters at different
    pipeline stages b/c variant calling for phasing happens before the integrative
    phasing (breakpointR/StrandPhaseR) may drop clusters. This cannot be resolved.
    --- unclear problem: it seems that WhatsHap produces phased VCFs for dropped clusters;
    presumably, there is still a function call somewhere that assumes consecutively numbered
    clusters.

    """
    import glob
    import sys
    import re

    glob_pattern = glob_path.format(**dict(wildcards))
    cluster_files = glob.glob(glob_pattern)
    if any('cluster99' in f for f in cluster_files):
        sys.stderr.write(f'\naux::check_cluster_file_completeness >>>\nWARNING: cluster99 file detected with glob pattern {glob_pattern}\n')
    cluster_files = set(f for f in cluster_files if 'cluster99' not in f)

    sample_name = sseq_reads.split('_')[0]
    if hasattr(wildcards, 'vc_reads'):
        readset = wildcards.vc_reads
    elif hasattr(wildcards, 'hap_reads'):
        readset = wildcards.hap_reads
    else:
        readset = None
    estimate_num_clusters, cluster_names = estimate_number_of_saarclusters(sample_name, wildcards.sseq_reads, return_names=True, readset=readset)
    num_clusters_slack = int(config.get('num_cluster_slack', 1))

    # first, check if there is one file per cluster name - estimate_number_of_saarclusters is accurate
    # after the integrative phasing stage of the pipeline
    produced_files = []
    for cn in cluster_names:
        find_name = re.compile(f'(\W){cn}(\W)')
        matched_file = [cluster_file for cluster_file in cluster_files if find_name.search(cluster_file) is not None]
        produced_files.extend(matched_file)

    if len(produced_files) == len(cluster_names):
        # we have a 1:1 match, this should be the best case
        cluster_files = produced_files
    elif not (estimate_num_clusters - num_clusters_slack) < len(cluster_files) < (estimate_num_clusters + num_clusters_slack):
        tmp = dict(wildcards)
        for cn in cluster_names:
            tmp[cluster_key] = cn
            cluster_file = source_path.format(**tmp)
            cluster_files.add(cluster_file)
    else:
        raise RuntimeError('Unconsidered situation')
    return sorted(cluster_files)


def check_cluster_no_variants(cluster_vcf_stats):

    no_variant_clusters = set()
    with open(cluster_vcf_stats, 'r') as table:
        for line in table:
            if 'HET-SNV_num' in line:
                record, number = line.strip().split()
                cluster_id = record.strip('_HET-SNV_num')
                if int(number) == 0:
                    no_variant_clusters.add(cluster_id)
    return no_variant_clusters


def check_cluster_vcf_is_empty(cluster_vcf_files):
    import pathlib as pl

    empty_clusters = set()
    for cluster_vcf_file in cluster_vcf_files:
        cluster_id = pl.Path(cluster_vcf_file).name.split('_')[0]
        is_empty = True
        with open(cluster_vcf_file, 'r') as vcf_file:
            for line in vcf_file:
                if line.startswith('#'):
                    continue
                is_empty = False
                assert line.strip()
                break
        if is_empty:
            empty_clusters.add(cluster_id)
    return empty_clusters


def collect_strandseq_alignments(wildcards, glob_collect=True, caller='snakemake'):
    """
    """
    import os
    import sys
    import glob

    func_name = '\nAGG::aux_utilities::collect_strandseq_alignments >> {}\n'

    if DEBUG:
        sys.stderr.write(func_name.format('wildcards ' + ' --- '.join('({}: {})'.format(k, v) for k, v in dict(wildcards).items())))

    source_path = os.path.join(
        'output',
        'alignments',
        'strandseq_to_reference',
        '{reference}',
        '{sseq_reads}',
        '{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.psort.mdup.sam.bam{ext}'
    )

    individual, project, platform_spec = wildcards.sseq_reads.split('_')[:3]
    platform, spec = platform_spec.split('-')

    glob_path = source_path.replace('{ext}', '*')
    glob_path = glob_path.replace('{lib_id}', '*')
    reduced_wildcards = dict(wildcards)
    reduced_wildcards['individual'] = individual
    reduced_wildcards['project'] = project
    reduced_wildcards['platform'] = platform
    reduced_wildcards['spec'] = spec

    pattern = glob_path.format(**reduced_wildcards)

    if DEBUG:
        sys.stderr.write(func_name.format('Glob collect w/ pattern: {}'.format(pattern)))
    
    sseq_libs, sseq_lib_ids = get_strandseq_library_info(wildcards.sseq_reads)

    if DEBUG:
        sys.stderr.write(func_name.format('Expecting N Strand-seq libraries: {}'.format(len(sseq_lib_ids))))

    bam_files = glob.glob(pattern)

    if len(bam_files) != len(sseq_lib_ids) * 2:  # factor 2: BAM plus BAM index
        if DEBUG:
            sys.stderr.write(func_name.format('Glob collect returned: {}'.format(len(bam_files))))
        
        bam_files = set(bam_files)

        for lib_id in sseq_lib_ids:
            reduced_wildcards['lib_id'] = lib_id
            reduced_wildcards['ext'] = ''
            bam_files.add(source_path.format(**reduced_wildcards))
            reduced_wildcards['ext'] = '.bai'
            bam_files.add(source_path.format(**reduced_wildcards))

    assert bam_files, 'collect_strandseq_alignments >> returned empty output: {}'.format(wildcards)

    if DEBUG:
        sys.stderr.write(func_name.format('Returning N BAM/BAM index files: {}'.format(len(bam_files))))

    return sorted(bam_files)


def get_sample_sex(sample):

    try:
        sample_desc = config['sample_description_{}'.format(sample)]
    except KeyError:
        raise ValueError('No sample description for {} in config - did you load the sample config YAML?'.format(sample))
    try:
        sample_sex = sample_desc['sex']
    except KeyError:
        sample_sex = 'unknown'
    return sample_sex


def load_saarclust_params_haploid(wildcards, input):
    return load_saarclust_params(wildcards, input, 'haploid')


def load_saarclust_params_squashed(wildcards, input):
    return load_saarclust_params(wildcards, input, 'squashed')


def load_saarclust_params(wildcards, input, use_case):

    key_map = dict((k, k.replace('_', '.')) for k in [
        'min_contig_size',
        'min_region_size',
        'bin_size',
        'step_size',
        'min_mapq',
        'allow_contig_cuts',
        'eval_ploidy'
    ])
    key_map['init_clusters'] = 'num.clusters'
    key_map['prob_threshold'] = 'prob.th'
    key_map['desired_clusters'] = 'desired.num.clusters'

    parameter_set = {
        'pairedReads': 'TRUE',
        'store.data.obj': 'TRUE',
        'reuse.data.obj': 'FALSE',
        'bin.method': '"dynamic"',
        'ord.method': '"greedy"',
        'assembly.fasta': '"{}"'.format(input.reference),
        'concat.fasta': None,
        'remove.always.WC': 'TRUE',
        'mask.regions': 'FALSE',
    }
    if use_case == 'haploid':
        parameter_set['concat.fasta'] = 'FALSE'
    elif use_case == 'squashed':
        parameter_set['concat.fasta'] = 'TRUE'
    else:
        raise ValueError('Unknown use case for SaaRclust parameter loading: {}'.format(use_case))
    
    for cfg_key, sc_key in key_map.items():
        if cfg_key == 'min_mapq':
            parameter_set[sc_key] = config.get(cfg_key, 0)
        else:
            parameter_set[sc_key] = config.get(cfg_key, None)

    individual = wildcards.sseq_reads.split('_')[0]
    sample_sex = get_sample_sex(individual)

    if sample_sex == 'male':
        parameter_set['desired.num.clusters'] = config.get('desired_clusters_male', config.get('desired_clusters', None))
        parameter_set['max.cluster.length.mbp'] = config.get('max_cluster_length_male_mbp', config.get('max_cluster_length_mbp', None))
    elif sample_sex == 'female':
        parameter_set['desired.num.clusters'] = config.get('desired_clusters_female', config.get('desired_clusters', None))
        parameter_set['max.cluster.length.mbp'] = config.get('max_cluster_length_female_mbp', config.get('max_cluster_length_mbp', None))
    else:
        parameter_set['desired.num.clusters'] = config.get('desired_clusters', None)

    use_case_to_rule = {
        'squashed': 'write_saarclust_config_file',
        'haploid': 'hac_write_saarclust_config_file'
    }
    
    non_default_params = config.get('sample_non_default_parameters', dict())
    if individual in non_default_params:
        sample_non_defaults = non_default_params[individual]
        use_non_defaults = True
        if 'use_only_in' in sample_non_defaults:
            try:
                sample_non_defaults = sample_non_defaults['use_only_in'][use_case_to_rule[use_case]]
            except KeyError:
                use_non_defaults = False

        if use_non_defaults:
            for cfg_key, sc_key in key_map.items():
                parameter_set[sc_key] = sample_non_defaults.get(cfg_key, parameter_set[sc_key])

    # drop all entries that were not specified in the pipeline config,
    # should default to whatever SaaRclust sets in this case
    parameter_set = dict((k, v) for k, v in parameter_set.items() if v is not None)

    # delete keys incompatible with earlier versions
    pipeline_version = int(config['git_commit_version'])
    if pipeline_version < 8:
        del parameter_set['desired.num.clusters']

    if pipeline_version < 9:
        del parameter_set['min.mapq']
    
    if pipeline_version < 14:
        del parameter_set['max.cluster.length.mbp']
        del parameter_set['allow.contig.cuts']
        del parameter_set['eval.ploidy']

    config_rows = ['[SaaRclust]'] + ['{} = {}'.format(k, parameter_set[k]) for k in sorted(parameter_set.keys())]

    saarclust_config = '\n'.join(config_rows) + '\n'
    return saarclust_config


def validate_checkpoint_output(check_output, expected_type='list_of_files'):
    """
    Because of github issues #55 and #142, this function exists to make
    sure the output returned from a checkpoint evaluation conforms to
    expectation (so far always list of file paths)

    :param check_output:
    :param expected_type:
    :return:
    """
    if expected_type == 'list_of_files':
        if isinstance(check_output, list):
            num_of_items = len(check_output)
            if num_of_items == 0:
                raise ValueError('Checkpoint evaluation resulted in empty list of files')
            for item in check_output:
                if os.path.isdir(item):
                    raise ValueError('Validating checkpoint output w/ {} items '
                                     '- encountered directory: {}'.format(num_of_items, item))
        elif isinstance(check_output, str):
            if os.path.isdir(check_output):
                # this is the reported case, handle in particular
                raise RuntimeError('Caught Snakemake error #55: checkpoint evaluation failed, '
                                   'returned single str / directory: {}'.format(check_output))
            else:
                raise ValueError('Single string returned from checkpoint evaluation, '
                                 'but is (not yet) a directory path: {}'.format(check_output))
        else:
            raise ValueError('Non-list type received for checkpoint '
                             'output validation: {} / {}'.format(type(check_output), str(check_output)))
    else:
        raise NotImplementedError('aux_utilities::validate_checkpoint_output called with '
                                  'unsupported expected type: {}'.format(expected_type))
    return


def load_preset_file(wildcards, input):
    """
    Save load for parameters from files that do not
    exist during a dry run, leading to a FileNotFoundError
    raised by Snakemake. Replaces construct

    params:
        preset = lambda wildcards, input: open(input.preset).read().strip()

    with

    params:
        preset = load_preset_file

    """
    if not hasattr(input, 'preset'):
        raise AttributeError('Input does not have a "preset" attribute: {}'.format(input))
    file_path = input.preset
    if not os.path.isfile(file_path):
        preset = 'PRESET-DRY-RUN'
    else:
        with open(file_path, 'r') as dump:
            preset = dump.read().strip()
            assert preset, 'Empty preset file: {}'.format(file_path)
    return preset


def load_fofn_file(input, prefix='', sep=' '):
    """
    Save load for list of filenames from a file that
    does not exist during a dry run, leading to a
    FileNotFoundError raised by Snakemake.

    Replaces construct

    params:
        preset = lambda wildcards, input: open(input.fofn).read().strip()

    with

    params:
        preset = lambda wildcards, input: load_fofn_file(input, prefix_string, separator_string)
    """
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


def load_seq_length_file(wildcards, input):
    """
    Follows same idea as two function above
    """
    if not hasattr(input, 'seq_info'):
        raise AttributeError('Input does not have "seq_info" attribute: {}'.format(input))
    file_path = input.seq_info
    if not os.path.isfile(file_path):
        seq_len = 'SEQLEN-DRY-RUN'
    else:
        seq_len = 0
        with open(file_path, 'r') as seq_info:
            for line in seq_info:
                if line.startswith('#'):
                    continue
                length = line.split()[1]
                try:
                    length = int(length)
                    seq_len += length
                except ValueError:
                    raise ValueError('Extracted seq. length is not an integer: {} / {}'.format(line.strip(), file_path))
    return seq_len


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
