
localrules: master_strandseq_dga_joint,
            write_haploid_split_reads_fofn

"""
Components:
vc_reads = FASTQ file used for variant calling relative to reference
hap_reads = FASTQ file to be used for haplotype reconstruction
sseq_reads = FASTQ file used for strand-seq phasing
"""
PATH_STRANDSEQ_DGA_JOINT = 'diploid_assembly/strandseq_joint/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sseq_reads}'
PATH_STRANDSEQ_DGA_JOINT_PROTECTED = PATH_STRANDSEQ_DGA_JOINT.replace('{', '{{').replace('}', '}}')


rule master_strandseq_dga_joint:
    input:
        []

def collect_haploid_split_reads(wildcards, glob_collect, file_ext):
    """
    :param wildcards:
    :param glob_collect:
    :param file_ext:
    :return:
    """
    import os

    read_format = {
        'fastq.gz': 'fastq',
        'pbn.bam': 'bam'
    }

    source_path = os.path.join('output',
                               PATH_STRANDSEQ_DGA_SPLIT,
                               'draft/haploid_' + read_format[file_ext],
                               '{hap_reads}.{hap}.{sequence}.' + file_ext)

    # pay attention to not merge the untagged reads twice here
    haplotypes = [wildcards.hap]
    if wildcards.hap.endswith('-un'):
        haplotypes = wildcards.hap.split('-')

    hap_read_splits = []

    if glob_collect:
        import glob
        for h in haplotypes:
            pattern = source_path.replace('{sequence}', '*')
            known_values = dict(wildcards)
            known_values['hap'] = h
            pattern = pattern.format(**known_values)
            hap_files = glob.glob(pattern)

            if not hap_files:
                raise RuntimeError('{}: collect_haploid_split_reads: no files collected with pattern {}'.format(h, pattern))
            hap_read_splits.extend(hap_files)

    else:
        reference_folder = os.path.join('output/reference_assembly/clustered', wildcards.sseq_reads)
        seq_output_dir = checkpoints.create_assembly_sequence_files.get(folder_path=reference_folder,
                                                                        reference=wildcards.reference).output[0]
        checkpoint_wildcards = glob_wildcards(os.path.join(seq_output_dir, '{sequence}.seq'))

        hap_read_splits = expand(
            source_path,
            var_caller=wildcards.var_caller,
            gq=wildcards.gq,
            qual=wildcards.qual,
            reference=wildcards.reference,
            vc_reads=wildcards.vc_reads,
            sseq_reads=wildcards.sseq_reads,
            hap_reads=wildcards.hap_reads,
            hap=haplotypes,  # note here: haplotype replacement
            sequence=checkpoint_wildcards.sequence
        )
    return hap_read_splits


def collect_haploid_split_reads_any(wildcards, glob_collect=False):
    """
    :param wildcards:
    :param glob_collect:
    :return:
    """
    if wildcards.subfolder == 'bam':
        collected_files = collect_haploid_split_reads(wildcards, glob_collect, 'pbn.bam')
    elif wildcards.subfolder == 'fastq':
        collected_files = collect_haploid_split_reads(wildcards, glob_collect, 'fastq.gz')
    else:
        raise ValueError('collect_haploid_split_reads_any: cannot process wildcards: {}'.format(wildcards))
    return sorted(collected_files)


rule write_haploid_split_reads_fofn:
    input:
        read_splits = collect_haploid_split_reads_any
    output:
        fofn = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_{subfolder}/{hap_reads}.{hap}.fofn',
    run:
        import os
        try:
            validate_checkpoint_output(input.read_splits)
            split_files = input.read_splits
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{}\n'.format(str(error)))
            split_files = collect_haploid_split_reads_any(wildcards, glob_collect=True)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(split_files):
                if not os.path.isfile(file_path):
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule concat_haploid_fastq:
    input:
        fofn = os.path.join('output', PATH_STRANDSEQ_DGA_JOINT, 'draft/haploid_fastq/{hap_reads}.{hap}.fofn')
    output:
        os.path.join('output', PATH_STRANDSEQ_DGA_JOINT, 'draft/haploid_fastq/{hap_reads}.{hap}.fastq.gz')
    conda:
         '../environment/conda/conda_shelltools.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        splits = lambda wildcards, input: load_fofn_file(input)
    shell:
         'cat {params.splits} > {output}'


rule concat_haploid_pbn_bam:
    input:
        fofn = os.path.join('output', PATH_STRANDSEQ_DGA_JOINT, 'draft/haploid_bam/{hap_reads}.{hap}.fofn')
    output:
        os.path.join('output', PATH_STRANDSEQ_DGA_JOINT, 'draft/haploid_bam/{hap_reads}.{hap}.pbn.bam')
    conda:
         '../environment/conda/conda_shelltools.yml'
    resources:
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    params:
        bam_parts = lambda wildcards, input: load_fofn_file(input, prefix=' -in ', sep=' -in ')
    shell:
        'bamtools merge {params.bam_parts} -out {output}'
