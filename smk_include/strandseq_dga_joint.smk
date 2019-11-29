
localrules: master_strandseq_dga_joint

"""
Components:
vc_reads = FASTQ file used for variant calling relative to reference
hap_reads = FASTQ file to be used for haplotype reconstruction
sts_reads = FASTQ file used for strand-seq phasing
"""
PATH_STRANDSEQ_DGA_JOINT = 'diploid_assembly/strandseq_joint/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}'


rule master_strandseq_dga_joint:
    input:


rule strandseq_dga_joint_haplo_tagging:
    input:
        vcf = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf.bgz',
        tbi = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{hap_reads}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{hap_reads}_map-to_{reference}.psort.sam.bam.bai',
        fasta = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta',
    output:
        bam = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haplotags/{hap_reads}.tagged.sam.bam',
        tags = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haplotags/{hap_reads}.tags.fq.tsv',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haplotags/{hap_reads}.tagging.fq.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haplotags/{hap_reads}.tagging.fq.rsrc'
    resources:
        runtime_hrs = 12,
        mem_total_mb = 4096,
        mem_per_cpu_mb = 4096
    shell:
        "whatshap --debug haplotag --output {output.bam} --reference {input.fasta} --output-haplotag-list {output.tags} {input.vcf} {input.bam} &> {log}"


rule strandseq_dga_joint_haplo_splitting:
    input:
        fastq = 'input/fastq/complete/{hap_reads}.fastq.gz',
        tags = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haplotags/{hap_reads}.tags.fq.tsv',
    output:
        h1 = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_fastq/{hap_reads}.h1.fastq.gz',
        h2 = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_fastq/{hap_reads}.h2.fastq.gz',
        un = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_fastq/{hap_reads}.un.fastq.gz',
        hist = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_fastq/{hap_reads}.rlen-hist.fq.tsv'
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_fastq/{hap_reads}.splitting.fq.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_fastq/{hap_reads}.splitting.fq.rsrc'
    resources:
        runtime_hrs = 12,
        mem_total_mb = 4096,
        mem_per_cpu_mb = 4096
    shell:
        "whatshap --debug split --pigz --output-h1 {output.h1} --output-h2 {output.h2} --output-untagged {output.un} --read-lengths-histogram {output.hist} {input.fastq} {input.tags} &> {log}"


rule strandseq_dga_joint_haplo_tagging_pacbio_native:
    input:
        vcf = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf.bgz',
        tbi = 'output/integrative_phasing/' + PATH_INTEGRATIVE_PHASING + '/{hap_reads}.wh-phased.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{hap_reads}_map-to_{reference}.psort.pbn.bam',
        bai = 'output/alignments/reads_to_reference/clustered/{sts_reads}/{hap_reads}_map-to_{reference}.psort.pbn.bam.bai',
        fasta = 'output/reference_assembly/clustered/{sts_reads}/{reference}.fasta',
    output:
        bam = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haplotags/{hap_reads}.tagged.pbn.bam',
        tags = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haplotags/{hap_reads}.tags.pbn.tsv',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haplotags/{hap_reads}.tagging.pbn.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haplotags/{hap_reads}.tagging.pbn.rsrc'
    resources:
        runtime_hrs = 12,
        mem_total_mb = 4096,
        mem_per_cpu_mb = 4096
    shell:
        "whatshap --debug haplotag --output {output.bam} --reference {input.fasta} --output-haplotag-list {output.tags} {input.vcf} {input.bam} &> {log}"


rule strandseq_dga_joint_haplo_splitting_pacbio_native:
    input:
        pbn_bam = 'input/bam/complete/{hap_reads}.pbn.bam',
        pbn_idx = 'input/bam/complete/{hap_reads}.pbn.bam.bai',
        tags = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haplotags/{hap_reads}.tags.pbn.tsv',
    output:
        h1 = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_bam/{hap_reads}.h1.pbn.bam',
        h2 = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_bam/{hap_reads}.h2.pbn.bam',
        un = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_bam/{hap_reads}.un.pbn.bam',
        hist = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_bam/{hap_reads}.rlen-hist.pbn.tsv'
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_bam/{hap_reads}.splitting.pbn.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_bam/{hap_reads}.splitting.pbn.rsrc',
    resources:
        runtime_hrs = 12,
        mem_total_mb = 4096,
        mem_per_cpu_mb = 4096
    shell:
        "whatshap --debug split --output-h1 {output.h1} --output-h2 {output.h2} --output-untagged {output.un} --read-lengths-histogram {output.hist} {input.pbn_bam} {input.tags} &> {log}"


rule strandseq_dga_joint_merge_tag_groups_pacbio_native:
    input:
        hap = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_bam/{hap_reads}.h{haplotype}.pbn.bam',
        un = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_bam/{hap_reads}.un.pbn.bam',
    output:
        'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_bam/{hap_reads}.h{haplotype}-un.pbn.bam',
    log:
        'log/output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_bam/{hap_reads}.h{haplotype}-un.pbn.mrg.log',
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_bam/{hap_reads}.h{haplotype}-un.pbn.mrg.rsrc',
    resources:
        runtime_hrs = 3
    wildcard_constraints:
        haplotype = '(1|2)'
    shell:
        'bamtools merge -in {input.hap} -in {input.un} -out {output} &> {log}'


rule strandseq_dga_joint_merge_tag_groups:
    input:
        hap = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_fastq/{hap_reads}.h{haplotype}.fastq.gz',
        un = 'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_fastq/{hap_reads}.un.fastq.gz',
    output:
        'output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_fastq/{hap_reads}.h{haplotype}-un.fastq.gz'
    benchmark:
        'run/output/' + PATH_STRANDSEQ_DGA_JOINT + '/draft/haploid_bam/{hap_reads}.h{haplotype}-un.fq.mrg.rsrc',
    resources:
        runtime_hrs = 3
    wildcard_constraints:
        haplotype = '(1|2)'
    shell:
        'cat {input.hap} {input.un} > {output}'
