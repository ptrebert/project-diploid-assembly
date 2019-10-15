
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'variant_calling.smk'

localrules: master_strandseq_dga

rule master_strandseq_dga:
    input:


rule install_rlib_breakpointr:
    input:
        'output/check_files/R_setup/saarclust.ok'
    output:
         'output/check_files/R_setup/breakpointr.ok'
    params:
        script_dir = config['script_dir']
    shell:
        'TAR=$(which tar) {params.script_dir}/install_breakpointr.R &> {output}'


rule install_rlib_strandphaser:
    input:
        'output/check_files/R_setup/breakpointr.ok'
    output:
         'output/check_files/R_setup/strandphaser.ok'
    params:
        script_dir = config['script_dir'],
        version = config['git_commit_strandphaser']
    shell:
        'TAR=$(which tar) {params.script_dir}/install_strandphaser.R {params.version} &> {output}'


rule write_breakpointr_config_file:
    input:
        setup_ok = 'output/check_files/R_setup/breakpointr.ok',
        reference = 'references/assemblies/{reference}.fasta'
    output:
        cfg = 'output/diploid_assembly/strandseq/config_files/{reference}/{sts_reads}/breakpointr.config'
    threads: config['num_cpu_high']
    run:
        config_rows = [
            '[General]',
            'numCPU = ' + str(threads),  # due to a bug in breakpointr, this value has to be repeated on the CLI
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


rule run_breakpointr:
    input:
        cfg = 'output/diploid_assembly/strandseq/config_files/{reference}/{sts_reads}/breakpointr.config',
        bams = collect_strandseq_alignments  # from module: prepare_custom_references - SaaRclust
    output:
        'output/diploid_assembly/strandseq/breakpointr/{reference}/{sts_reads}/{reference}.WCregions.txt',
        'output/diploid_assembly/strandseq/breakpointr/{reference}/{sts_reads}/breakpointR.config',
        'output/diploid_assembly/strandseq/breakpointr/{reference}/{sts_reads}/README.txt',
        directory('output/diploid_assembly/strandseq/breakpointr/{reference}/{sts_reads}/breakpoints'),
        directory('output/diploid_assembly/strandseq/breakpointr/{reference}/{sts_reads}/browserfiles'),
        directory('output/diploid_assembly/strandseq/breakpointr/{reference}/{sts_reads}/data'),
        directory('output/diploid_assembly/strandseq/breakpointr/{reference}/{sts_reads}/plots')
    log:
        'log/output/diploid_assembly/strandseq/breakpointr/{reference}/{sts_reads}/breakpointr.log'
    benchmark:
        'run/output/diploid_assembly/strandseq/breakpointr/{reference}/{sts_reads}/breakpointr.rsrc'
    params:
        input_dir = lambda wildcards, input: os.path.dirname(input.bams[0]),
        output_dir = lambda wildcards, output: os.path.dirname(output[0]),
        script_dir = config['script_dir']
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 128,
        mem_total_mb = config['num_cpu_high'] * 128
    shell:
        '{params.script_dir}/run_breakpointr.R {params.input_dir} {input.cfg} {params.output_dir} {threads} {output} &> {log}'


rule write_strandphaser_config_file:
    input:
        setup_ok = 'output/check_files/R_setup/strandphaser.ok',
        reference = 'references/assemblies/{reference}.fasta'
    output:
        cfg = 'output/diploid_assembly/strandseq/config_files/{reference}/{sts_reads}/strandphaser.config'
    threads: config['num_cpu_high']
    run:
        config_rows = [
            '[General]',
            'numCPU = ' + str(threads),
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


rule run_strandphaser:
    input:
        cfg = 'output/diploid_assembly/strandseq/config_files/{reference}/{sts_reads}/strandphaser.config',
        bams = collect_strandseq_alignments,  # from module: prepare_custom_references - SaaRclust
        wc_regions = 'output/diploid_assembly/strandseq/breakpointr/{reference}/{sts_reads}/{reference}.WCregions.txt',
        variant_calls = 'output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{vc_reads}.final.vcf'
    output:
        browser = directory('output/strandseq_phased_variants/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/browserFiles'),
        data = directory('output/strandseq_phased_variants/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/data'),
        phased = directory('output/strandseq_phased_variants/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/Phased'),
        maps = directory('output/strandseq_phased_variants/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/SingleCellHaps'),
        vcf_dir = directory('output/strandseq_phased_variants/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/VCFfiles'),
        cfg = 'output/strandseq_phased_variants/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/StrandPhaseR.config',
        vcf = 'output/strandseq_phased_variants/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/final/{sts_reads}.phased.vcf'
    log:
        stp = 'log/output/strandseq_phased_variants/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}.phased.log',
        bcf = 'log/output/strandseq_phased_variants/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}.concat.log',
    benchmark:
        'run/output/strandseq_phased_variants/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}.phased.rsrc'
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 256,
        mem_total_mb = config['num_cpu_high'] * 256
    wildcard_constraints:
        gq = '[0-9]+',
        dp = '[0-9]+',
        reference = '[\w\-]+',
        vc_reads = '[\w\-]+',
        sts_reads = '[\w\-]+'
    params:
        input_dir = lambda wildcards, input: os.path.dirname(input.bams[0]),
        output_dir = lambda wildcards, output: os.path.dirname(output.cfg),
        individual = lambda wildcards: wildcards.sts_reads.split('_')[0],
        script_dir = config['script_dir']
    shell:
        '{params.script_dir}/run_strandphaser.R {params.input_dir} {input.cfg} ' \
            ' {input.variant_calls} {input.wc_regions} ' \
            ' {params.output_dir} {params.individual} &> {log.stp}' \
            ' && ' \
            ' bcftools concat --output {output.vcf} --output-type v `ls {output.vcf_dir}/*phased.vcf` ' \
            ' &> {log.bcf}'


rule strandseq_dga_phase_variants:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{vc_reads}.final.vcf.bgz',
        tbi = 'output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{vc_reads}.final.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.sam.bam.bai',
        fasta = 'references/assemblies/{reference}.fasta',
        seq_info = 'references/assemblies/{reference}/sequences/{sequence}.seq',
        sts_phased = 'output/strandseq_phased_variants/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/final/{sts_reads}.phased.vcf'
    output:
        vcf = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.phased.vcf'
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.phased.log'
    benchmark:
        'run/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.phased.rsrc'
    wildcard_constraints:
        var_caller = '(freebayes|longshot)',
        gq = '[0-9]+',
        dp = '[0-9]+',
        reference = '[\w\-]+',
        vc_reads = '[\w\-]+',
        sts_reads = '[\w\-]+'
    shell:
        'whatshap --debug phase --chromosome {wildcards.sequence} --reference {input.fasta} ' \
            ' {input.vcf} {input.bam} {input.sts_phased} 2> {log} ' \
            ' | egrep "^(#|{wildcards.sequence}\s)" > {output}'


def sdga_collect_sequence_phased_vcf_files(wildcards):
    """
    """
    seq_output_dir = checkpoints.create_assembly_sequence_files.get(**wildcards).output[0]

    checkpoint_wildcards = glob_wildcards(
        os.path.join(seq_output_dir, '{sequence}.seq')
        )

    vcf_files = expand(
        'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/split_by_seq/{hap_reads}.{sequence}.phased.vcf',
        var_caller=wildcards.var_caller,
        reference=wildcards.reference,
        gq=wildcards.gq,
        dp=wildcards.dp,
        vc_reads=wildcards.vc_reads,
        sts_reads=wildcards.sts_reads,
        hap_reads=wildcards.hap_reads,
        sequence=checkpoint_wildcards.sequence
        )
    return sorted(vcf_files)


rule strandseq_dga_merge_sequence_phased_vcf_files:
    input:
        vcf_files = sdga_collect_sequence_phased_vcf_files
    output:
        'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.phased.vcf'
    wildcard_constraints:
        gq = '[0-9]+',
        dp = '[0-9]+',
        reference = '[\w\-]+',
        vc_reads = '[\w\-]+',
        sts_reads = '[\w\-]+',
        hap_reads = '[\w\-]+'
    shell:
        'bcftools concat --output {output} --output-type v {input.vcf_files}'


rule strandseq_dga_haplo_tagging:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        vcf = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.phased.vcf.bgz',
        tbi = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.phased.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.sam.bam.bai',
        fasta = 'references/assemblies/{reference}.fasta'
    output:
        bam = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.tagged.sam.bam',
        tags = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.tags.tsv',
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.tagging.log',
    benchmark:
        'run/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.tagging.rsrc',
    conda:
        config['conda_env_whsplit']
    shell:
        "whatshap --debug haplotag --output {output.bam} --reference {input.fasta} --output-haplotag-list {output.tags} {input.vcf} {input.bam} &> {log}"


rule strandseq_dga_haplo_splitting:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        fastq = 'input/fastq/complete/{hap_reads}.fastq.gz',
        tags = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.tags.tsv',
    output:
        h1 = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.h1.fastq.gz',
        h2 = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.h2.fastq.gz',
        un = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.un.fastq.gz',
        hist = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.rlen-hist.tsv'
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.splitting.log',
    benchmark:
        'run/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.splitting.rsrc',
    conda:
        config['conda_env_whsplit']
    wildcard_constraints:
        gq = '[0-9]+',
        dp = '[0-9]+',
        reference = '[\w\-]+',
        vc_reads = '[\w\-]+',
        sts_reads = '[\w\-]+',
        hap_reads = '[\w\-]+'
    shell:
        "whatshap --debug split --pigz --output-h1 {output.h1} --output-h2 {output.h2} --output-untagged {output.un} --read-lengths-histogram {output.hist} {input.fastq} {input.tags} &> {log}"


rule strandseq_dga_haplo_tagging_pacbio_native:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = PacBio native BAM file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        vcf = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.phased.vcf.bgz',
        tbi = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.phased.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.pbn.bam',
        bai = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.pbn.bam.bai',
        fasta = 'references/assemblies/{reference}.fasta'
    output:
        bam = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.tagged.pbn.bam',
        tags = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.tags.pbn.tsv',
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.tagging.pbn.log',
    benchmark:
        'run/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.tagging.pbn.rsrc',
    conda:
        config['conda_env_whsplit']
    shell:
        "whatshap --debug haplotag --output {output.bam} --reference {input.fasta} --output-haplotag-list {output.tags} {input.vcf} {input.bam} &> {log}"


rule strandseq_dga_haplo_splitting_pacbio_native:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = PacBio native BAM file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        pbn_bam = 'input/bam/complete/{hap_reads}.pbn.bam',
        pbn_idx = 'input/bam/complete/{hap_reads}.pbn.bam.bai',
        tags = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.tags.pbn.tsv',
    output:
        h1 = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.h1.pbn.bam',
        h2 = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.h2.pbn.bam',
        un = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.un.pbn.bam',
        hist = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.rlen-hist.tsv'
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.splitting.pbn.log',
    benchmark:
        'run/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.splitting.pbn.rsrc',
    conda:
        config['conda_env_whsplit']
    wildcard_constraints:
        gq = '[0-9]+',
        dp = '[0-9]+',
        reference = '[\w\-]+',
        vc_reads = '[\w\-]+',
        sts_reads = '[\w\-]+',
        hap_reads = '[\w\-]+'
    shell:
        "whatshap --debug split --pigz --output-h1 {output.h1} --output-h2 {output.h2} --output-untagged {output.un} --read-lengths-histogram {output.hist} {input.pbn_bam} {input.tags} &> {log}"


rule strandseq_dga_merge_tag_groups_pacbio_native:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = PacBio native BAM file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        hap = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.h{hap}.pbn.bam',
        un = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.un.pbn.bam',
    output:
        'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.h{hap}-un.pbn.bam',
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.h{hap}-un.pbn.mrg.log',
    wildcard_constraints:
        hap = '(1|2)'
    shell:
        'bamtools merge -in {input.hap} -in {input.un} -out {output} &> {log}'


rule strandseq_dga_merge_tag_groups:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        hap = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.h{haplotype}.fastq.gz',
        un = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.un.fastq.gz',
    output:
        'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.h{haplotype}-un.fastq.gz',
    wildcard_constraints:
        haplotype = '(1|2)'
    shell:
        'cat {input.hap} {input.un} > {output}'


rule strandseq_dga_assemble_haplotypes_wtdbg_layout:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        fastq = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}.{hap}.fastq.gz',
    output:
        layout = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/layout/wtdbg2/{hap_reads}.{hap}.ctg.lay.gz',
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}-wtdbg.{hap}.layout.log',
    benchmark:
        'run/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{hap_reads}-wtdbg.{hap}.layout.rsrc',
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 6144,
        mem_total_mb = 409600
    params:
        param_preset = lambda wildcards: config['wtdbg2_presets'][wildcards.hap_reads.rsplit('_', 1)[0]],
        out_prefix = lambda wildcards, output: output.layout.rsplit('.', 3)[0]
    shell:
        'wtdbg2 -x {params.param_preset} -i {input.fastq} -g3g -t {threads} -o {params.out_prefix} &> {log}'


rule strandseq_dga_assemble_haplotypes_wtdbg_consensus:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sts_reads = FASTQ file used for strand-seq phasing
    """
    input:
        layout = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/layout/wtdbg2/{hap_reads}.{hap}.ctg.lay.gz',
    output:
        'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}-wtdbg.{hap}.fasta',
    wildcard_constraints:
        hap = '[h12\-un]+'
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}-wtdbg.{hap}.log',
    benchmark:
        'run/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}-wtdbg.{hap}.rsrc',
    threads: config['num_cpu_high']
    resources:
        mem_per_cpu_mb = 384,
        mem_total_mb = 12288
    wildcard_constraints:
        gq = '[0-9]+',
        dp = '[0-9]+',
        reference = '[\w\-]+',
        vc_reads = '[\w\-]+',
        sts_reads = '[\w\-]+',
        hap_reads = '[\w\-]+',
        hap = '[h12un\-]+'
    shell:
        'wtpoa-cns -t {threads} -i {input.layout} -o {output} &> {log}'
