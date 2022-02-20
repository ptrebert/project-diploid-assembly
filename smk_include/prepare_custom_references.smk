
localrules: master_prepare_custom_references,
            write_saarclust_config_file,


rule master_prepare_custom_references:
    input:
        []


def collect_strandseq_merge_files(wildcards, glob_collect=True, caller='snakemake'):
    """
    Replacing checkpoints with regular rules creates inherent problems for Strand-seq
    data because the potential library QC selection makes it impossible to know
    which libraries will eventually be used. However, this should not be a problem here
    because the dual fraction Strand-seq samples are never channeled through automatic QC.
    """
    import glob
    import fnmatch

    source_path = 'output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/aln/{individual}_{project}_{platform}-{spec}_{lib_id}_{run_id}.filt.sam.bam'

    individual = wildcards.individual
    platform = wildcards.platform
    project = wildcards.project
    lib_id = wildcards.lib_id

    assert individual in wildcards.reference and individual in wildcards.sseq_reads, \
        'Wrong reference / sseq_reads match: {} / {}'.format(wildcards.reference, wildcards.sseq_reads)

    pattern = source_path.replace('{run_id}', '*')
    pattern = pattern.replace('{spec}', '*')
    pattern = pattern.format(**dict(wildcards))
    bam_files = glob.glob(pattern)

    if not bam_files or len(bam_files) != 2:  # note here: no BAM index files considered
        # use info from module "scrape_data_sources" to determine rule input
        sseq_libs, sseq_lib_ids = get_strandseq_library_info(wildcards.sseq_reads)
        selected_libs = fnmatch.filter(sseq_libs, '*{}*'.format(wildcards.lib_id))
        if len(selected_libs) != 2:
            raise RuntimeError('Cannot identify Strand-seq libs for merge: {} / {}'.format(wildcards, selected_libs))
        new_source_path = 'output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/aln/{sseq_lib}.filt.sam.bam'
        bam_files = []
        for lib in selected_libs:
            tmp = dict(wildcards)
            tmp['sseq_lib'] = lib
            bam_files.append(new_source_path.format(**tmp))

    assert len(bam_files) == 2, 'Unexpected number of BAM files for merge: {}'.format(bam_files)

    return sorted(bam_files)


rule merge_mono_dinucleotide_fraction:
    """
    In the process of removing all checkpoints from the pipeline,
    move the input collection directly into this rule
    Previous input:
        fofn = rules.write_strandseq_merge_fofn.output.fofn
    """
    input:
        bams = collect_strandseq_merge_files
    output:
        temp('output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.sam.bam')
    log:
        'log/output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.log'
    benchmark:
        'rsrc/output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/mrg/{individual}_{project}_{platform}-{spec}_{lib_id}.mrg.rsrc'
    wildcard_constraints:
        sseq_reads = CONSTRAINT_STRANDSEQ_DIFRACTION_SAMPLES
    conda:
        '../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_low']
    shell:
        'samtools merge -@ {threads} -O BAM {output} {input.bams} &> {log}'


rule link_strandseq_monofraction_samples:
    """
    Switch to copying here because sym linking seems to cause
    trouble (not recognized) on certain file systems
    """
    input:
        'output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/aln/{library_id}.filt.sam.bam'
    output:
        'output/alignments/strandseq_to_reference/{reference}/{sseq_reads}/temp/mrg/{library_id}.mrg.sam.bam'
    wildcard_constraints:
        sseq_reads = CONSTRAINT_STRANDSEQ_MONOFRACTION_SAMPLES
    shell:
        'cp {input} {output}'


rule samtools_position_sort_strandseq_reads:
    """
    Since Strand-seq alignments are small, make dedicated
    samtools sort rule with lower resource requirements
    """
    input:
        '{folder_path}/temp/mrg/{sts_library}.mrg.sam.bam'
    output:
        temp('{folder_path}/temp/sort/{sts_library}.mrg.psort.sam.bam')
    conda:
        '../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sts_library = '[A-Za-z0-9\-_]+'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_total_mb = lambda wildcards, attempt: config['num_cpu_low'] * 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'samtools sort -m {resources.mem_per_cpu_mb}M --threads {threads} -o {output} {input}'


rule mark_duplicate_reads_strandseq:
    input:
        rules.samtools_position_sort_strandseq_reads.output[0]
    output:
        '{folder_path}/{sts_library}.mrg.psort.mdup.sam.bam'
    log:
        'log/{folder_path}/{sts_library}.mrg.psort.mdup.log'
    benchmark:
        'rsrc/{folder_path}/{sts_library}.mrg.psort.mdup.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        sts_library = '[A-Za-z0-9\-_]+'
    threads: config['num_cpu_low']
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_total_mb = lambda wildcards, attempt: config['num_cpu_low'] * 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    shell:
        'sambamba markdup -t {threads} --overflow-list-size 600000 {input} {output} &> {log}'


rule write_saarclust_config_file:
    """
    As long as aggregate-style input functions downstream of
    checkpoints are problematic in a cluster environment,
    avoid complications by making the config writing local
    (aggregate should work), and explicitly write a file
    containing just the input folder for StrandPhaseR
    """
    input:
        setup_ok = rules.install_rlib_saarclust.output.check,
        reference = 'output/reference_assembly/non-hap-res/{reference}.fasta',
        ref_idx = 'output/reference_assembly/non-hap-res/{reference}.fasta.fai',
        strandseq_reads = 'input/fastq/{sseq_reads}.fofn',
        bam = collect_strandseq_alignments  # from module aux_utilities
    output:
        cfg = 'output/reference_assembly/clustered/temp/saarclust/config/{reference}/{sseq_reads}/saarclust.config',
        input_dir = 'output/reference_assembly/clustered/temp/saarclust/config/{reference}/{sseq_reads}/saarclust.input'
    params:
        saarclust = lambda wildcards, input: load_saarclust_params(wildcards, input, 'squashed')
    run:
        import os

        outfolder = os.path.dirname(list(input.bam)[0])
        assert os.path.isdir(outfolder), 'write_saarclust_config_file: alignment output folder does not exist: {}'.format(outfolder)

        with open(output.cfg, 'w') as dump:
            _ = dump.write(params.saarclust)

        with open(output.input_dir, 'w') as dump:
            _ = dump.write(outfolder + '\n')


rule run_saarclust_assembly_clustering:
    """
    Converted from checkpoint to regular rule to get rid of
    checkpoint-related problems.
    This one is quite problematic b/c SaaRclust does not produce
    a fix number of clusters.
    """
    input:
        cfg = rules.write_saarclust_config_file.output.cfg,
        fofn = rules.write_saarclust_config_file.output.input_dir
    output:
        ctg_report = 'output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/clustered_assembly/ctgReport_2e+05bp_dynamic.tsv',
        dir_data = directory('output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/data'),
        dir_plots = directory('output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/plots'),
        cfg = 'output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/SaaRclust.config',
    log:
        'log/output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/saarclust.log'
    benchmark:
        'rsrc/output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/saarclust.rsrc'
    conda:
        '../environment/conda/conda_rscript.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 + 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 + 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    params:
        script_exec = lambda wildcards: find_script_path('run_saarclust.R'),
        out_folder = lambda wildcards, output: os.path.dirname(output.cfg),
        in_folder = lambda wildcards, input: load_fofn_file(input),
        version = config['git_commit_saarclust']
    shell:
        '{params.script_exec} {input.cfg} {params.in_folder} {params.out_folder} {params.version} &> {log} '


def collect_clustered_fasta_sequences(wildcards, glob_collect=False, caller='snakemake'):
    """
    """
    raise RuntimeError('collect_clustered_fasta_sequences: function is deprecated')
    import sys
    debug = bool(config.get('show_debug_messages', False))
    warn = bool(config.get('show_warnings', False))
    func_name = '\nCALL::{}\nchk::agg::collect_clustered_fasta_sequences: {{}}\n'.format(caller)

    if debug:
        sys.stderr.write(func_name.format('wildcards ' + str(wildcards)))

    source_path = 'output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/clustered_assembly/{sequence}.fasta'

    if glob_collect:
        if debug:
            sys.stderr.write(func_name.format('called w/ glob collect'))
        import glob
        pattern = source_path.replace('{sequence}', '*')
        pattern = pattern.format(**dict(wildcards))
        fasta_files = glob.glob(pattern)

        if not fasta_files:
            if debug:
                sys.stderr.write(func_name.format('glob collect failed'))
            raise RuntimeError('collect_clustered_fasta_sequences: no files collected with pattern {}'.format(pattern))

    else:
        raise RuntimeError('Illegal function call: Snakemake checkpoints must not be used.')
        if debug:
            sys.stderr.write(func_name.format('called w/ chk::get'))
        from snakemake.exceptions import IncompleteCheckpointException as ICE

        strandseq_reads = wildcards.sseq_reads
        nhr_assembly = wildcards.reference

        try:
            # this output folder is the /clustered_assembly subfolder
            seq_output_dir = checkpoints.run_saarclust_assembly_clustering.get(reference=nhr_assembly, sseq_reads=strandseq_reads).output.dir_fasta
            checkpoint_wildcards = glob_wildcards(os.path.join(seq_output_dir, '{sequence}.fasta'))

            fasta_files = expand(
                source_path,
                reference=nhr_assembly,
                sseq_reads=strandseq_reads,
                sequence=checkpoint_wildcards.sequence
            )
        except ICE as ice:
            if debug:
                sys.stderr.write(func_name.format('chk::get raised SMK::ICE'))
            try:
                fasta_files = collect_clustered_fasta_sequences(wildcards, glob_collect=True, caller='debug-internal')
            except RuntimeError:
                if debug:
                    sys.stderr.write(func_name.format('glob collect failed - re-raising SMK::ICE'))
                raise ice
            else:
                if debug or warn:
                    sys.stderr.write(func_name.format('Snakemake error: glob collect success, but SMK::ICE raised (#216)'))
        else:
            if debug:
                sys.stderr.write(func_name.format('chk::get did not raise ICE - checkpoint passed'))

    return sorted(fasta_files)


# DEPRECATED
# rule write_reference_fasta_clusters_fofn:
#     """
#     Local rule with minimal overhead - properly collect checkpoint output
#     """
#     input:
#         fasta = collect_clustered_fasta_sequences
#     output:
#         fofn = 'output/reference_assembly/clustered/temp/saarclust/{sseq_reads}/{reference}.clusters.fofn'
#     run:
#         try:
#             validate_checkpoint_output(input.fasta)
#             fasta_files = input.fasta
#         except (RuntimeError, ValueError) as error:
#             import sys
#             sys.stderr.write('\n{}\n'.format(str(error)))
#             fasta_files = collect_clustered_fasta_sequences(wildcards, glob_collect=True, caller='write_reference_fasta_clusters_fofn')

#         with open(output.fofn, 'w') as dump:
#             for file_path in sorted(fasta_files):
#                 if not os.path.isfile(file_path):
#                     import sys
#                     sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
#                 _ = dump.write(file_path + '\n')


rule run_gba_rescue:
    """
    For contigs that cannot be clustered by Strand-seq,
    bring them back into the assembly via guilt-by-association
    (check assembly graph)

    2022-02-20
    for the hg002/acro experiments, it was necessary to add additional
    sequences into the NHR assembly before clustering. The script below
    was adapted to also keep these "unknown" sequences in the stats output,
    but only so if the new parameter "--ignore-non-assembly-sequence" is set;
    otherwise, the script will raise and report unknown sequence in the assembly.
    """
    input:
        ctg_report = 'output/reference_assembly/clustered/temp/saarclust/results/{hap_reads}_nhr-{assembler}/{sseq_reads}/clustered_assembly/ctgReport_2e+05bp_dynamic.tsv',
        raw_unitigs = 'output/reference_assembly/non-hap-res/layout/{assembler}/{hap_reads}/{hap_reads}.r_utg.noseq.gfa',
        p_contigs = 'output/reference_assembly/non-hap-res/layout/{assembler}/{hap_reads}/{hap_reads}.p_ctg.noseq.gfa',
    output:
        table = 'output/statistics/clustering/{sseq_reads}/{hap_reads}_nhr-{assembler}.cluster-info.tsv',
        unassigned = 'output/statistics/clustering/{sseq_reads}/{hap_reads}_nhr-{assembler}.unassigned-unitigs.tsv',
    conda:
        '../environment/conda/conda_pyscript.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    params:
        script_exec = lambda wildcards: find_script_path('connect_unitig_contig.py'),
    shell:
        '{params.script_exec} -u {input.raw_unitigs} -c {input.p_contigs} -s {input.ctg_report} '
            '-o {output.table} --unassigned-unitigs {output.unassigned} '


def load_contig_fasta_sequences(file_path, is_oriented):
    """
    Load and cache contig sequences; if the input is not
    the raw NHR assembly (= contig sequences are oriented
    using Strand-seq information), the contig information is extracted
    from position 1 after splitting the FASTA header line at "_"
    """
    cache = dict()
    this_seq = ''
    this_ctg = ''
    with open(file_path, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                if this_ctg:
                    cache[this_ctg] = this_seq
                this_ctg = line.strip()[1:]
                if is_oriented:
                    # position 0 is cluster ID
                    this_ctg = this_ctg.split('_')[1]
                this_seq = ''
            else:
                this_seq += line.strip()
    cache[this_ctg] = this_seq
    return cache


rule write_clustered_split_fasta:
    input:
        table = 'output/statistics/clustering/{sseq_reads}/{hap_reads}_nhr-{assembler}.cluster-info.tsv',
        nhr_fasta = 'output/reference_assembly/non-hap-res/{hap_reads}_nhr-{assembler}.fasta',
    output:
        fasta = 'output/reference_assembly/clustered/{{sseq_reads}}/{{hap_reads}}_scV{}-{{assembler}}.fasta-split.fa'.format(config['git_commit_version'])
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    run:
        import pandas as pd
        import collections as col
        import io
        nhr_seq = load_contig_fasta_sequences(input.nhr_fasta, is_oriented=False)
        df = pd.read_csv(input.table, sep='\t', header=0)
        # input_cluster_sizes
        init_cluster_sizes = df.drop_duplicates(['cluster_id', 'contig'], inplace=False).groupby('cluster_id')['contig_length'].sum().to_dict()
        REVCOMP_TABLE = get_revcomp_translation_table()

        sample = wildcards.hap_reads.split('_')[0]
        sample_sex = get_sample_sex(sample)[0].upper()
        garbage = io.StringIO()
        dump_garbage = False
        dumped_cluster_sizes = col.Counter()
        seen = set()
        with open(output.fasta, 'w') as dump:
            for (cluster_id, contig, scl_key), unitigs in df.groupby(['cluster_id', 'contig', 'scl_key']):
                if contig in seen:
                    raise RuntimeError(f'Duplicate contig: {input.table} / {contig}')
                output_name = f'>{cluster_id}_{contig}_{sample}-{sample_sex}_SCL-{scl_key}'
                output_seq = nhr_seq[contig]
                if 'dR' in scl_key:
                    # directionality of contig is revcomp
                    output_seq = output_seq[::-1].translate(REVCOMP_TABLE)
                dumped_cluster_sizes[cluster_id] += len(output_seq)
                if cluster_id == 'cluster99':
                    _ = garbage.write(f'{output_name}\n{output_seq}\n')
                    dump_garbage = True
                else:
                    _ = dump.write(f'{output_name}\n{output_seq}\n')
                seen.add(contig)
        
        size_errors = []
        for cluster, size_expect in init_cluster_sizes.items():
            if cluster == 'cluster99':
                continue
            size_observed = dumped_cluster_sizes[cluster]
            if size_expect != size_observed:
                size_errors.append((cluster, size_expect, size_observed))

        if size_errors:
            raise RuntimeError(f'fasta-split cluster size errors: {size_errors}')

        if dump_garbage:
            cluster99_out = output.fasta.replace('fasta-split', 'cluster99')
            with open(cluster99_out, 'w') as dump:
                _ = dump.write(garbage.getvalue())
    # END OF RUN BLOCK


def dump_cluster_size_error(sseq_reads, contig, seq_len, max_seq_len, ignore_size_error):

    if seq_len > max_seq_len:
        with open(f'ERROR_cluster-size_{sseq_reads}.err', 'w'):
            pass
        if not ignore_size_error:
            raise ValueError(f'Squashed assembly cluster too large ({contig} - {seq_len} / max {max_seq_len}): {sseq_reads}')
    return


rule write_clustered_concat_fasta:
    input:
        table = 'output/statistics/clustering/{sseq_reads}/{hap_reads}_nhr-{assembler}.cluster-info.tsv',
        split_fasta = 'output/reference_assembly/clustered/{{sseq_reads}}/{{hap_reads}}_scV{}-{{assembler}}.fasta-split.fa'.format(config['git_commit_version'])
    output:
        fasta = 'output/reference_assembly/clustered/{{sseq_reads}}/{{hap_reads}}_scV{}-{{assembler}}.fasta'.format(config['git_commit_version']),
        cluster_ids = 'output/reference_assembly/clustered/{{sseq_reads}}/{{hap_reads}}_scV{}-{{assembler}}.cluster-ids.txt'.format(config['git_commit_version'])
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    run:
        import pandas as pd
        import collections as col
        import io
        # by using the contigs of the fasta-split assembly, it is ensured
        # that the directionality information is propagated to the concatenated
        # assembly FASTA
        ctg_seq = load_contig_fasta_sequences(input.split_fasta, is_oriented=True)
        df = pd.read_csv(input.table, sep='\t', header=0)
        init_cluster_sizes = df.drop_duplicates(['cluster_id', 'contig'], inplace=False).groupby('cluster_id')['contig_length'].sum().to_dict()

        MAX_SEQ_LENGTH = 268435456
        IGNORE_SIZE_ERROR = bool(config.get('ignore_cluster_size_error', False))

        spacer = 'N' * 100
        last_cluster = ''
        last_seq = ''
        last_length = 0
        dumped_clusters = []
        dumped_cluster_sizes = col.Counter()
        with open(output.fasta, 'w') as dump:
            for (cluster_id, contig), unitigs in df.groupby(['cluster_id', 'contig']):
                if cluster_id == 'cluster99':
                    continue
                if cluster_id != last_cluster:
                    if last_seq:
                        dump_cluster_size_error(
                            wildcards.sseq_reads, last_cluster,
                            last_length, MAX_SEQ_LENGTH,
                            IGNORE_SIZE_ERROR
                        )
                        _ = dump.write(f'>{last_cluster}\n{last_seq}\n')
                        dumped_clusters.append(last_cluster)

                    last_cluster = cluster_id
                    last_seq = ctg_seq[contig]
                    last_length = len(last_seq)
                    dumped_cluster_sizes[last_cluster] += last_length
                else:
                    this_seq = ctg_seq[contig]
                    last_seq += spacer + this_seq
                    last_length += len(this_seq)
                    dumped_cluster_sizes[last_cluster] += len(this_seq)
            if last_seq:
                _ = dump.write(f'>{last_cluster}\n{last_seq}\n')
                dumped_clusters.append(last_cluster)

        size_errors = []
        for cluster, size_expect in init_cluster_sizes.items():
            if cluster == 'cluster99':
                continue
            size_observed = dumped_cluster_sizes[cluster]
            if size_expect != size_observed:
                size_errors.append((cluster, size_expect, size_observed))

        if size_errors:
            raise RuntimeError(f'fasta-concat cluster size errors: {size_errors}')

        with open(output.cluster_ids, 'w') as dump:
            _ = dump.write('\n'.join(dumped_clusters) + '\n')
    # END OF RUN BLOCK


# rule write_reference_fasta_clusters_fofn:
#     """
#     Collect output of SaaRclust-ering
#     NB: reference here is "nhr" assembly

#     The "clusters.fofn" output is processed internally - if available - in the function
#     aux_utilities::estimate_number_of_saarclusters

#     """
#     input:
#         fasta_dir = 'output/reference_assembly/clustered/temp/saarclust/results/{reference}/{sseq_reads}/clustered_assembly',
#         fasta_ref = 'output/reference_assembly/non-hap-res/{reference}.fasta',
#     output:
#         fofn = 'output/reference_assembly/clustered/temp/saarclust/{sseq_reads}/{reference}.clusters.fofn'
#     run:
#         import os
#         import sys

#         show_warnings = bool(config.get('show_warnings', False))

#         if not os.path.isdir(input.fasta_dir):
#             raise RuntimeError('SaaRclust output folder does not exist: {}'.format(input.fasta_dir))

#         cluster_files = []
#         for file_name in os.listdir(input.fasta_dir):
#             if not file_name.endswith('.fasta'):
#                 continue
#             if 'cluster99' in file_name:
#                 if show_warnings:
#                     sys.stderr.write('SaaRclust / cluster99 detected for sample {} / {}\n'.format(wildcards.reference, wildcards.sseq_reads))
#                 continue
#             file_path = os.path.join(input.fasta_dir, file_name)
#             cluster_files.append(file_path)
        
#         if not cluster_files:
#             raise RuntimeError('No FASTA cluster files detected in SaaRclust output directory: {}'.format(input.fasta_dir))
        
#         with open(output.fofn, 'w') as fofn:
#             _ = fofn.write('\n'.join(sorted(cluster_files)))


# rule check_max_cluster_size:
#     """
#     This is a sanity check to avoid that downstream contig alignment tasks fail
#     due to squashed assembly clusters being too large to be processed (restriction of SAM/BAM)
#     see:
#     https://github.com/lh3/minimap2/issues/440#issuecomment-508052956
#     """
#     input:
#         'output/reference_assembly/clustered/temp/saarclust/{sseq_reads}/{reference}.clusters.fofn'
#     output:
#         'output/reference_assembly/clustered/temp/saarclust/{sseq_reads}/{reference}.clusters.size.ok'
#     run:
#         max_seq_len = 268435456
#         cache = dict()
#         ignore_size_error = bool(config.get('ignore_cluster_size_error', False))
#         with open(input[0], 'r') as fasta_list:
#             for fasta_file in fasta_list:
#                 seq_len = 0
#                 cluster_name = ''
#                 with open(fasta_file.strip(), 'r') as fasta:
#                     for line in fasta:
#                         if line.startswith('>'):
#                             cluster_name = line.strip().strip('>')
#                             continue
#                         seq_len += len(line.strip())
#                 if seq_len > max_seq_len:
#                     # create an untracked error output file for simpler
#                     # summary of what samples failed and why
#                     with open(f'ERROR_cluster-size_{wildcards.sseq_reads}.err', 'w'):
#                         pass
#                     if not ignore_size_error:
#                         raise ValueError('Squashed assembly cluster too large ({} / max {}): {}'.format(seq_len, max_seq_len, fasta_file))
#                 cache[cluster_name] = seq_len

#         with open(output[0], 'w') as check_ok:
#             for c in sorted(cache.keys()):
#                 cluster_size = cache[c]
#                 ratio = round(cluster_size / max_seq_len * 100, 2)
#                 _ = check_ok.write('{}\t{}\t{}\n'.format(c, cluster_size, ratio))


# rule merge_reference_fasta_clusters:
#     input:
#         fofn = 'output/reference_assembly/clustered/temp/saarclust/{sseq_reads}/{hap_reads}_nhr-{assembler}.clusters.fofn',
#         size_ok = 'output/reference_assembly/clustered/temp/saarclust/{sseq_reads}/{hap_reads}_nhr-{assembler}.clusters.size.ok'
#     output:
#         expand('output/reference_assembly/clustered/{{sseq_reads}}/{{hap_reads}}_scV{version}-{{assembler}}.fasta',
#                 version=config['git_commit_version'])
#     conda:
#         '../environment/conda/conda_shelltools.yml'
#     params:
#         fasta_clusters = lambda wildcards, input: load_fofn_file(input)
#     shell:
#         'cat {params.fasta_clusters} > {output}'
