
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'variant_calling.smk'

localrules: master_canonical_dga

rule master_canonical_dga:
    input:


rule canonical_dga_phase_variants:
    """
    vc_readset = FASTQ file used for variant calling relative to reference
    hap_readset = FASTQ file to be used for haplotype reconstruction
    """
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{vc_reads}.final.vcf.bgz',
        tbi = 'output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{vc_reads}.final.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.sam.bam.bai',
        fasta = 'references/assemblies/{reference}.fasta'
    output:
        vcf = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.phased.vcf'
    log:
        'log/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.phased.log'
    benchmark:
        'run/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.phased.rsrc'
    shell:
        'whatshap --debug phase --reference {input.fasta} --output {output.vcf} {input.vcf} {input.bam} &> {log}'


rule canonical_dga_haplo_tagging:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    """
    input:
        vcf = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.phased.vcf.bgz',
        tbi = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.phased.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/{hap_reads}_map-to_{reference}.psort.sam.bam.bai',
        fasta = 'references/assemblies/{reference}.fasta'
    output:
        bam = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.tagged.sam.bam',
        tags = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.tags.tsv',
    log:
        'log/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.tagging.log',
    benchmark:
        'run/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.tagging.rsrc',
    conda:
        "../environment/conda/wh_split.yml"
    shell:
        "whatshap --debug haplotag --output {output.bam} --reference {input.fasta} --output-haplotag-list {output.tags} {input.vcf} {input.bam} &> {log}"


rule canonical_dga_haplo_splitting:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    """
    input:
        fastq = 'input/fastq/complete/{hap_reads}.fastq.gz',
        tags = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.tags.tsv',
    output:
        h1 = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.h1.fastq.gz',
        h2 = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.h2.fastq.gz',
        un = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.un.fastq.gz',
        hist = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.rlen-hist.tsv'
    log:
        'log/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.splitting.log',
    benchmark:
        'run/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.splitting.rsrc',
    conda:
        "../environment/conda/wh_split.yml"
    shell:
        "whatshap --debug split --pigz --output-h1 {output.h1} --output-h2 {output.h2} --output-untagged {output.un} --read-lengths-histogram {output.hist} {input.fastq} {input.tags} &> {log}"


rule canonical_dga_merge_tag_groups:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    """
    input:
        hap = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.h{hap}.fastq.gz',
        un = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.un.fastq.gz',
    output:
        'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.h{hap}-un.fastq.gz',
    wildcard_constraints:
        hap = '(1|2)'
    shell:
        'cat {input.hap} {input.un} > {output}'


rule canonical_dga_assemble_haplotypes_layout:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    """
    input:
        fastq = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.{hap}.fastq.gz',
    output:
        layout = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/layout/{hap_reads}.{hap}.ctg.lay.gz',
    wildcard_constraints:
        haplotype = '[h12\-un]+'
    log:
        'log/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.{hap}.layout.log',
    benchmark:
        'run/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.{hap}.layout.rsrc',
    threads: 64
    params:
        param_preset = lambda wildcards: config['wtdbg2_presets'][wildcards.hap_reads.rsplit('_', 1)[0]],
        out_prefix = lambda wildcards, output: output.layout.rsplit('.', 3)[0]
    shell:
        'wtdbg2 -x {params.param_preset} -i {input.fastq} -g3g -t {threads} -o {params.out_prefix} &> {log}'


rule canonical_dga_assemble_haplotypes_consensus:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    """
    input:
        layout = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/layout/{hap_reads}.{hap}.ctg.lay.gz',
    output:
        'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/consensus/{hap_reads}.{hap}.fasta',
    wildcard_constraints:
        hap = '[h12\-un]+'
    log:
        'log/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/consensus/{hap_reads}.{hap}.log',
    benchmark:
        'run/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/consensus/{hap_reads}.{hap}.rsrc',
    threads: 64
    shell:
        'wtpoa-cns -t {threads} -i {input.layout} -o {output} &> {log}'