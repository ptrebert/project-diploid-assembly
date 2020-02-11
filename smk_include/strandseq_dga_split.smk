
localrules: master_strandseq_dga_split


"""
Components:
vc_reads = FASTQ file used for variant calling relative to reference
hap_reads = FASTQ file to be used for haplotype reconstruction
sts_reads = FASTQ file used for strand-seq phasing
"""
PATH_STRANDSEQ_DGA_SPLIT = 'diploid_assembly/strandseq_split/{var_caller}_QUAL{qual}_GQ{gq}/{reference}/{vc_reads}/{sts_reads}'
PATH_STRANDSEQ_DGA_SPLIT_PROTECTED = PATH_STRANDSEQ_DGA_SPLIT.replace('{', '{{').replace('}', '}}')


rule master_strandseq_dga_split:
    input:
        []


rule strandseq_dga_split_haplo_tagging:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    fq_hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        vcf = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{fq_hap_reads}.wh-phased.vcf.bgz',
        tbi = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{fq_hap_reads}.wh-phased.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{fq_hap_reads}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{fq_hap_reads}_map-to_{reference}.psort.sam.bam.bai',
        fasta = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta',
        seq_info = 'output/reference_assembly/clustered/{sts_reads}/{reference}/sequences/{sequence}.seq',
    output:
        bam = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{fq_hap_reads}.{sequence}.tagged.sam.bam',
        tags = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{fq_hap_reads}.{sequence}.tags.fq.tsv',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{fq_hap_reads}.{sequence}.fq.tag.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{fq_hap_reads}.{sequence}.fq.tag.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        fq_hap_reads = CONSTRAINT_ALL_FASTQ_INPUT_SAMPLES + '_[0-9]+',
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    shell:
        'whatshap --debug haplotag --regions {wildcards.sequence} --output {output.bam} '
            '--reference {input.fasta} --output-haplotag-list {output.tags} '
            '{input.vcf} {input.bam} &> {log}'


rule strandseq_dga_split_haplo_tagging_pacbio_native:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    pbn_hap_reads = PacBio native BAM file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        vcf = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{pbn_hap_reads}.wh-phased.vcf.bgz',
        tbi = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{pbn_hap_reads}.wh-phased.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{pbn_hap_reads}_map-to_{reference}.psort.pbn.bam',
        bai = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{pbn_hap_reads}_map-to_{reference}.psort.pbn.bam.bai',
        fasta = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta',
        seq_info = 'output/reference_assembly/clustered/{sts_reads}/{reference}/sequences/{sequence}.seq',
    output:
        bam = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{pbn_hap_reads}.{sequence}.tagged.pbn.bam',
        tags = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{pbn_hap_reads}.{sequence}.tags.pbn.tsv',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{pbn_hap_reads}.{sequence}.pbn.tag.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{pbn_hap_reads}.{sequence}.pbn.tag.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        pbn_hap_reads = CONSTRAINT_ALL_PBN_INPUT_SAMPLES + '_[0-9]+',
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: 4 * attempt
    shell:
        'whatshap --debug haplotag --regions {wildcards.sequence} --output {output.bam} '
            '--reference {input.fasta} --output-haplotag-list {output.tags} '
            '{input.vcf} {input.bam} &> {log}'


rule strandseq_dga_split_haplo_splitting:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    fq_hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        fastq = 'input/fastq/complete/{fq_hap_reads}.fastq.gz',
        tags = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{fq_hap_reads}.{sequence}.tags.fq.tsv',
    output:
        h1 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{fq_hap_reads}.h1.{sequence}.fastq.gz',
        h2 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{fq_hap_reads}.h2.{sequence}.fastq.gz',
        un = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{fq_hap_reads}.un.{sequence}.fastq.gz',
        hist = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{fq_hap_reads}.rlen-hist.{sequence}.fq.tsv'
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{fq_hap_reads}.{sequence}.fq.split.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{fq_hap_reads}.{sequence}.fq.split.rsrc',
    conda:
        '../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        fq_hap_reads = CONSTRAINT_ALL_FASTQ_INPUT_SAMPLES + '_[0-9]+',
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    shell:
        'whatshap --debug split --discard-unknown-reads --pigz --only-largest-block '
            '--output-h1 {output.h1} --output-h2 {output.h2} --output-untagged {output.un} '
            '--read-lengths-histogram {output.hist} {input.fastq} {input.tags} &> {log}'


rule strandseq_dga_split_haplo_splitting_pacbio_native:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    pbn_hap_reads = PacBio native BAM file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        pbn_bam = 'input/bam/complete/{pbn_hap_reads}.pbn.bam',
        pbn_idx = 'input/bam/complete/{pbn_hap_reads}.pbn.bam.bai',
        tags = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haplotags/{pbn_hap_reads}.{sequence}.tags.pbn.tsv',
    output:
        h1 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pbn_hap_reads}.h1.{sequence}.pbn.bam',
        h2 = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pbn_hap_reads}.h2.{sequence}.pbn.bam',
        un = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pbn_hap_reads}.un.{sequence}.pbn.bam',
        hist = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pbn_hap_reads}.rlen-hist.{sequence}.pbn.tsv'
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pbn_hap_reads}.{sequence}.pbn.split.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pbn_hap_reads}.{sequence}.pbn.split.rsrc',
    conda:
        '../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        pbn_hap_reads = CONSTRAINT_ALL_PBN_INPUT_SAMPLES + '_[0-9]+',
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    shell:
        'whatshap --debug split --discard-unknown-reads --only-largest-block '
            '--output-h1 {output.h1} --output-h2 {output.h2} --output-untagged {output.un} '
            '--read-lengths-histogram {output.hist} {input.pbn_bam} {input.tags} &> {log}'


rule strandseq_dga_split_merge_tag_groups:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    fq_hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        hap = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{fq_hap_reads}.h{haplotype}.{sequence}.fastq.gz',
        un = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{fq_hap_reads}.un.{sequence}.fastq.gz',
    output:
        'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{fq_hap_reads}.h{haplotype}-un.{sequence}.fastq.gz'
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_fastq/{fq_hap_reads}.h{haplotype}-un.{sequence}.fq.mrg.rsrc'
    conda:
        '../environment/conda/conda_shelltools.yml'
    wildcard_constraints:
        haplotype = '(1|2)',
        fq_hap_reads = CONSTRAINT_ALL_FASTQ_INPUT_SAMPLES + '_[0-9]+',
    shell:
        'cat {input.hap} {input.un} > {output}'


rule strandseq_dga_split_merge_tag_groups_pacbio_native:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    pbn_hap_reads = PacBio native BAM file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        hap = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pbn_hap_reads}.h{haplotype}.{sequence}.pbn.bam',
        un = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pbn_hap_reads}.un.{sequence}.pbn.bam',
    output:
        'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pbn_hap_reads}.h{haplotype}-un.{sequence}.pbn.bam',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pbn_hap_reads}.h{haplotype}-un.{sequence}.pbn.mrg.log'
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_bam/{pbn_hap_reads}.h{haplotype}-un.{sequence}.pbn.mrg.rsrc'
    conda:
        '../environment/conda/conda_biotools.yml'
    wildcard_constraints:
        haplotype = '(1|2)',
        pbn_hap_reads = CONSTRAINT_ALL_PBN_INPUT_SAMPLES + '_[0-9]+',
    resources:
        mem_total_mb = lambda wildcards, attempt: 512 + 512 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 512 + 512 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'bamtools merge -in {input.hap} -in {input.un} -out {output} &> {log}'


# The following merge step is primarily for convenience, i.e. for pushing
# the assembled sequence as single FASTA file into the evaluation part
# of the pipeline (e.g., QUAST-LG)

def collect_assembled_sequence_files(wildcards, glob_collect=False):
    """
    """
    import os

    source_path = os.path.join('output',
                               PATH_STRANDSEQ_DGA_SPLIT,
                               'draft/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.fasta')

    if glob_collect:
        import glob
        pattern = source_path.replace('{sequence}', '*')
        pattern = pattern.format(**dict(wildcards))
        fasta_files = glob.glob(pattern)

        if not fasta_files:
            raise RuntimeError('collect_assembled_sequence_files: no files collected with pattern {}'.format(pattern))

    else:
        reference_folder = os.path.join('output/reference_assembly/clustered', wildcards.sts_reads)
        seq_output_dir = checkpoints.create_assembly_sequence_files.get(folder_path=reference_folder,
                                                                        reference=wildcards.reference).output[0]
        checkpoint_wildcards = glob_wildcards(os.path.join(seq_output_dir, '{sequence}.seq'))

        fasta_files = expand(
            source_path,
            var_caller=wildcards.var_caller,
            gq=wildcards.gq,
            qual=wildcards.qual,
            reference=wildcards.reference,
            vc_reads=wildcards.vc_reads,
            sts_reads=wildcards.sts_reads,
            hap_reads=wildcards.hap_reads,
            assembler=wildcards.assembler,
            hap=wildcards.hap,
            sequence=checkpoint_wildcards.sequence
        )
    return fasta_files


rule write_assembled_fasta_clusters_fofn:
    input:
        cluster_fastas = collect_assembled_sequence_files
    output:
        fofn = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-{assembler}.{hap}.fofn',
    run:
        import os
        try:
            validate_checkpoint_output(input.cluster_fastas)
            fasta_files = input.cluster_fastas
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{}\n'.format(str(error)))
            fasta_files = collect_assembled_sequence_files(wildcards, glob_collect=True)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(fasta_files):
                if not os.path.isfile(file_path):
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule strandseq_dga_split_merge_assembled_cluster_fastas:
    """
    """
    input:
        fofn = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-{assembler}.{hap}.fofn'
    output:
         'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/draft/haploid_assembly/{hap_reads}-{assembler}.{hap}.fasta'
    params:
        cluster_fastas = lambda wildcards, input: load_fofn_file(input)
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 * attempt
    run:
        import sys
        import io
        import re
        import collections
        fasta_header = '^>[A-Za-z_0-9]+$'

        header_counter = collections.Counter()

        with open(output[0], 'w') as merged_fasta:
            pass

        for fasta_file in params.cluster_fastas.split():
            buffer = io.StringIO()
            _, seq_id, _ = fasta_file.rsplit('.', 2)
            with open(fasta_file, 'r') as single_fasta:
                for line in single_fasta:
                    if line.startswith('>'):
                        line = line.replace('>', '>{}_'.format(seq_id))
                        # make sure there is no bogus info in the header
                        line = line.split()[0].strip()
                        if re.match(fasta_header, line) is None:
                            # implies other characters in header
                            sys.stderr.write('\nERROR: malformed FASTA header detected during merge '
                                             'of assembled cluster FASTA files - affected file {} / '
                                             'header {}\n'.format(fasta_file, line))
                            try:
                                os.unlink(output[0])
                            except (OSError, IOError):
                                # can't do anything about this
                                pass
                            raise ValueError('Malformed FASTA header in merge: {}'.format(line))
                        header_counter[line] += 1
                        line += '\n'
                    buffer.write(line)

            with open(output[0], 'a') as merged_fasta:
                _ = merged_fasta.write(buffer.getvalue())

        for header, count in header_counter.most_common():
            if count <= 1:
                break
            sys.stderr.write('\nERROR: duplicate header in assembled FASTA merge: {} / {}\n'.format(header, count))
            if os.path.isfile(output[0]):
                # note: deleting output leads to failed job
                os.unlink(output[0])


def collect_polished_contigs(wildcards, glob_collect=False):
    """
    :param wildcards:
    :param glob_collect:
    :return:
    """
    import os

    source_path = os.path.join('output',
                               PATH_STRANDSEQ_DGA_SPLIT,
                               'polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{sequence}.{pol_pass}.fasta')

    if glob_collect:
        import glob
        pattern = source_path.replace('{sequence}', '*')
        pattern = pattern.format(**dict(wildcards))
        fasta_files = glob.glob(pattern)

        if not fasta_files:
            raise RuntimeError('collect_polished_contigs: no files collected with pattern {}'.format(pattern))

    else:
        reference_folder = os.path.join('output/reference_assembly/clustered', wildcards.sts_reads)
        seq_output_dir = checkpoints.create_assembly_sequence_files.get(folder_path=reference_folder,
                                                                        reference=wildcards.reference).output[0]

        checkpoint_wildcards = glob_wildcards(
            os.path.join(seq_output_dir, '{sequence}.seq')
            )

        fasta_files = expand(
            source_path,
            var_caller=wildcards.var_caller,
            qual=wildcards.qual,
            gq=wildcards.gq,
            reference=wildcards.reference,
            vc_reads=wildcards.vc_reads,
            sts_reads=wildcards.sts_reads,
            pol_reads=wildcards.pol_reads,
            hap_reads=wildcards.hap_reads,
            assembler=wildcards.assembler,
            hap=wildcards.hap,
            sequence=checkpoint_wildcards.sequence,
            pol_pass=wildcards.pol_pass
        )

    return fasta_files


rule write_polished_contigs_fofn:
    input:
        contigs = collect_polished_contigs
    output:
        fofn = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{pol_pass}.fofn'
    run:
        import os

        try:
            validate_checkpoint_output(input.contigs)
            fasta_files = input.contigs
        except (RuntimeError, ValueError) as error:
            import sys
            sys.stderr.write('\n{]\n'.format(str(error)))
            fasta_files = collect_polished_contigs(wildcards, glob_collect=True)

        with open(output.fofn, 'w') as dump:
            for file_path in sorted(fasta_files):
                if not os.path.isfile(file_path):
                    import sys
                    sys.stderr.write('\nWARNING: File missing, may not be created yet - please check: {}\n'.format(file_path))
                _ = dump.write(file_path + '\n')


rule merge_polished_contigs:
    input:
        fofn = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{pol_pass}.fofn'
    output:
        fofn = 'output/' + PATH_STRANDSEQ_DGA_SPLIT + '/polishing/{pol_reads}/haploid_assembly/{hap_reads}-{assembler}.{hap}.{pol_pass}.fasta'
    params:
        polished_fastas = lambda wildcards, input: load_fofn_file(input)
    resources:
        mem_total_mb = lambda wildcards, attempt: 1024 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 1024 * attempt
    run:
        import sys
        import io
        import re
        import collections
        fasta_header = '^>[A-Za-z_0-9]+$'

        header_counter = collections.Counter()

        with open(output[0], 'w') as merged_fasta:
            pass

        for fasta_file in params.polished_fastas.split():
            buffer = io.StringIO()
            _, seq_id, _ , _ = fasta_file.rsplit('.', 3)
            with open(fasta_file, 'r') as single_fasta:
                for line in single_fasta:
                    if line.startswith('>'):
                        if seq_id not in line:
                            line = line.replace('>', '>{}_'.format(seq_id))
                        # some awesome workaround for silly names coming from other tools...
                        line = line.replace('|arrow', '')
                        # make sure there is no other bogus info in the header
                        line = line.split()[0].strip()
                        if re.match(fasta_header, line) is None:
                            # implies other characters in header
                            sys.stderr.write('\nERROR: malformed FASTA header detected during merge '
                                             'of assembled cluster FASTA files - affected file {} / '
                                             'header {}\n'.format(fasta_file, line))
                            try:
                                os.unlink(output[0])
                            except (OSError, IOError):
                                # can't do anything about this
                                pass
                            raise ValueError('Malformed FASTA header in merge: {}'.format(line))
                        header_counter[line] += 1
                        line += '\n'
                    buffer.write(line)
            with open(output[0], 'a') as merged_fasta:
                _ = merged_fasta.write(buffer.getvalue())

        for header, count in header_counter.most_common():
            if count <= 1:
                break
            sys.stderr.write('\nERROR: duplicate header in polished FASTA merge: {} / {}\n'.format(header, count))
            if os.path.isfile(output[0]):
                # note: deleting output leads to failed job
                os.unlink(output[0])
