
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'integrative_phasing.smk'
include: 'variant_calling.smk'

localrules: master_strandseq_dga_split, \
            strandseq_dga_split_merge_readsets, \
            strandseq_dga_split_merge_sequences


"""
Components:
vc_reads = FASTQ file used for variant calling relative to reference
hap_reads = FASTQ file to be used for haplotype reconstruction
sts_reads = FASTQ file used for strand-seq phasing
"""
PATH_STRANDSEQ_DGA_SPLIT = 'diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}'


rule master_strandseq_dga_split:
    input:


rule strandseq_dga_split_haplo_tagging:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        vcf = 'output/integrative_phasing/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.phased.vcf.bgz',
        tbi = 'output/integrative_phasing/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.phased.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.sam.bam.bai',
        fasta = 'references/assemblies/{reference}.fasta',
        seq_info = 'references/assemblies/{reference}/sequences/{sequence}.seq',
    output:
        bam = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.tagged.sam.bam',
        tags = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.tags.fq.tsv',
    log:
        'log/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.tagging.fq.log',
    benchmark:
        'run/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.tagging.fq.rsrc',
    shell:
        "whatshap --debug haplotag --regions {sequence} --output {output.bam} --reference {input.fasta} --output-haplotag-list {output.tags} {input.vcf} {input.bam} &> {log}"


rule strandseq_dga_split_haplo_tagging_pacbio_native:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = PacBio native BAM file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        vcf = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.phased.vcf.bgz',
        tbi = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.phased.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.pbn.bam',
        bai = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.pbn.bam.bai',
        fasta = 'references/assemblies/{reference}.fasta',
        seq_info = 'references/assemblies/{reference}/sequences/{sequence}.seq',
    output:
        bam = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.tagged.pbn.bam',
        tags = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.tags.pbn.tsv',
    log:
        'log/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.tagging.pbn.log',
    benchmark:
        'run/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.tagging.pbn.rsrc',
    shell:
        "whatshap --debug haplotag --regions {sequence} --output {output.bam} --reference {input.fasta} --output-haplotag-list {output.tags} {input.vcf} {input.bam} &> {log}"


rule strandseq_dga_split_haplo_splitting:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        fastq = 'input/fastq/complete/{hap_reads}.fastq.gz',
        tags = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.tags.fq.tsv',
    output:
        h1 = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.h1.{sequence}.fastq.gz',
        h2 = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.h2.{sequence}.fastq.gz',
        un = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.un.{sequence}.fastq.gz',
        hist = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.rlen-hist.{sequence}.fq.tsv'
    log:
        'log/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.splitting.fq.log',
    benchmark:
        'run/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.splitting.fq.rsrc',
    resources:
        mem_per_cpu_mb = 8192,
        mem_total_mb = 8192
    shell:
        "whatshap --debug split --discard-unknown-reads --pigz --output-h1 {output.h1} --output-h2 {output.h2} --output-untagged {output.un} --read-lengths-histogram {output.hist} {input.fastq} {input.tags} &> {log}"


rule strandseq_dga_split_haplo_splitting_pacbio_native:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = PacBio native BAM file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        pbn_bam = 'input/bam/complete/{hap_reads}.pbn.bam',
        pbn_idx = 'input/bam/complete/{hap_reads}.pbn.bam.bai',
        tags = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.tags.pbn.tsv',
    output:
        h1 = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.h1.{sequence}.pbn.bam',
        h2 = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.h2.{sequence}.pbn.bam',
        un = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.un.{sequence}.pbn.bam',
        hist = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.rlen-hist.{sequence}.pbn.tsv'
    log:
        'log/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.splitting.pbn.log',
    benchmark:
        'run/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.splitting.pbn.rsrc',
    resources:
        mem_per_cpu_mb = 8192,
        mem_total_mb = 8192
    shell:
        "whatshap --debug split --discard-unknown-reads --output-h1 {output.h1} --output-h2 {output.h2} --output-untagged {output.un} --read-lengths-histogram {output.hist} {input.pbn_bam} {input.tags} &> {log}"


rule strandseq_dga_split_merge_tag_groups:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        hap = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.h{haplotype}.{sequence}.fastq.gz',
        un = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.un.{sequence}.fastq.gz',
    output:
        'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.h{haplotype}-un.{sequence}.fastq.gz',
    wildcard_constraints:
        haplotype = '(1|2)'
    shell:
        'cat {input.hap} {input.un} > {output}'


rule strandseq_dga_split_merge_tag_groups_pacbio_native:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = PacBio native BAM file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        hap = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.h{haplotype}.{sequence}.pbn.bam',
        un = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.un.{sequence}.pbn.bam',
    output:
        'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.h{haplotype}-un.{sequence}.pbn.bam',
    log:
        'log/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.h{haplotype}-un.{sequence}.pbn.mrg.log',
    wildcard_constraints:
        haplotype = '(1|2)'
    shell:
        'bamtools merge -in {input.hap} -in {input.un} -out {output} &> {log}'


rule strandseq_dga_split_assemble_haplotypes_wtdbg_layout:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        fastq = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{hap}.{sequence}.fastq.gz',
        seq_info = 'references/assemblies/{reference}/sequences/{sequence}.seq',
    output:
        layout = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/layout/wtdbg2/{sequence}/{hap_reads}.{hap}.ctg.lay.gz',
    log:
        'log/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}-wtdbg.{hap}.{sequence}.layout.log',
    benchmark:
        'run/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}-wtdbg.{hap}.{sequence}.layout.rsrc',
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 409600
    params:
        param_preset = lambda wildcards: config['wtdbg2_presets'][wildcards.hap_reads.rsplit('_', 1)[0]],
        out_prefix = lambda wildcards, output: output.layout.rsplit('.', 3)[0],
        seq_len = lambda wildcards, input: open(input.seq_info).readline().split()[1]
    shell:
        'wtdbg2 -x {params.param_preset} -i {input.fastq} -g {params.seq_len} -t {threads} -o {params.out_prefix} &> {log}'


rule strandseq_dga_split_assemble_haplotypes_wtdbg_consensus:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        layout = 'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/layout/wtdbg2/{sequence}/{hap_reads}.{hap}.ctg.lay.gz',
    output:
        'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/consensus/{hap_reads}-wtdbg.{hap}.{sequence}.fasta',
    log:
        'log/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/consensus/{hap_reads}-wtdbg.{hap}.{sequence}.log',
    benchmark:
        'run/output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/consensus/{hap_reads}-wtdbg.{hap}.{sequence}.rsrc',
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 384,
        mem_total_mb = 12288
    shell:
        'wtpoa-cns -t {threads} -i {input.layout} -o {output} &> {log}'


def collect_sequence_readsets(wildcards):
    """
    """
    seq_output_dir = checkpoints.create_assembly_sequence_files.get(reference=wildcards.reference).output[0]

    checkpoint_wildcards = glob_wildcards(seq_output_dir, '{sequence}.seq')

    seq_files = expand('output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{hap}.{sequence}.fastq.gz',
                        var_caller=wildcards.var_caller,
                        gq=wildcards.gq,
                        dp=wildcards.dp,
                        reference=wildcards.reference,
                        vc_reads=wildcards.vc_reads,
                        sts_reads=wildcards.sts_reads,
                        hap_reads=wildcards.hap_reads,
                        assembler=wildcards.assembler,
                        hap=wildcards.hap,
                        sequence=checkpoint_wildcards.sequence)
    return sorted(seq_files)


rule strandseq_dga_split_merge_readsets:
    """
    """
    input:
        sequence_readsets = collect_sequence_readsets
    output:
        'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.fastq.gz',
    shell:
        'cat {input} > {output}'

def collect_assembled_sequence_files(wildcards):
    """
    """
    seq_output_dir = checkpoints.create_assembly_sequence_files.get(reference=wildcards.reference).output[0]

    checkpoint_wildcards = glob_wildcards(seq_output_dir, '{sequence}.seq')

    seq_files = expand('output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/consensus/{hap_reads}-{assembler}.{hap}.{sequence}.fasta',
                        var_caller=wildcards.var_caller,
                        gq=wildcards.gq,
                        dp=wildcards.dp,
                        reference=wildcards.reference,
                        vc_reads=wildcards.vc_reads,
                        sts_reads=wildcards.sts_reads,
                        hap_reads=wildcards.hap_reads,
                        assembler=wildcards.assembler,
                        hap=wildcards.hap,
                        sequence=checkpoint_wildcards.sequence)
    return sorted(seq_files)


rule strandseq_dga_split_merge_sequences:
    """
    """
    input:
        assembled_sequences = collect_assembled_sequence_files
    output:
        'output/diploid_assembly/strandseq_split/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}-{assembler}.{hap}.fasta',
    shell:
        'cat {input} > {output}'