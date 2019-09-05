
include: 'aux_utilities.smk'
include: 'preprocess_input.smk'
include: 'preprocess_references.smk'
include: 'prepare_custom_references.smk'
include: 'variant_calling.smk'

localrules: master_canonical_dga

rule master_canonical_dga:
    input:
        # runs on CORDELIA
        expand('output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/consensus/{readset}.{hap}.fasta',
                var_caller=['longshot'],
                gq=config['filter_vcf_gq'],
                dp=config['filter_vcf_dp'],
                reference=['HG00733_hgsvc_pbsq2-ccs_sqa-100kb'],
                readset=['HG00733_hgsvc_pbsq2-ccs_1000'],
                hap=['h1', 'h2', 'h1-un', 'h2-un'])


rule canonical_dga_phase_variants:
    input:
        vcf = 'output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{readset}.final.vcf.bgz',
        tbi = 'output/variant_calls/{var_caller}/{reference}/final_GQ{gq}_DP{dp}/{readset}.final.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/{readset}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/{readset}_map-to_{reference}.psort.sam.bam.bai',
        fasta = 'references/assemblies/{reference}.fasta'
    output:
        vcf = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.phased.vcf'
    log:
        'log/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.phased.log'
    benchmark:
        'run/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.phased.rsrc'
    shell:
        'whatshap --debug phase --reference {input.fasta} --output {output.vcf} {input.vcf} {input.bam} &> {log}'


rule canonical_dga_haplo_tagging:
    input:
        vcf = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.phased.vcf.bgz',
        tbi = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.phased.vcf.bgz.tbi',
        bam = 'output/alignments/reads_to_reference/{readset}_map-to_{reference}.psort.sam.bam',
        bai = 'output/alignments/reads_to_reference/{readset}_map-to_{reference}.psort.sam.bam.bai',
        fasta = 'references/assemblies/{reference}.fasta'
    output:
        bam = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.tagged.sam.bam',
        tags = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.tags.tsv',
    log:
        'log/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.tagging.log',
    benchmark:
        'run/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.tagging.rsrc',
    conda:
        "../environment/conda/wh_split.yml"
    shell:
        "whatshap --debug haplotag --output {output.bam} --reference {input.fasta} --output-haplotag-list {output.tags} {input.vcf} {input.bam} &> {log}"


rule canonical_dga_haplo_splitting:
    input:
        fastq = 'input/fastq/complete/{readset}.fastq.gz',
        tags = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.tags.tsv',
    output:
        h1 = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.h1.fastq.gz',
        h2 = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.h2.fastq.gz',
        un = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.un.fastq.gz',
        hist = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.rlen-hist.tsv'
    log:
        'log/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.splitting.log',
    benchmark:
        'run/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.splitting.rsrc',
    conda:
        "../environment/conda/wh_split.yml"
    shell:
        "whatshap --debug split --pigz --output-h1 {output.h1} --output-h2 {output.h2} --output-untagged {output.un} --read-lengths-histogram {output.hist} {input.fastq} {input.tags} &> {log}"


rule canonical_dga_merge_tag_groups:
    input:
        hap = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.h{hap}.fastq.gz',
        un = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.un.fastq.gz',
    output:
        'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.h{hap}-un.fastq.gz',
    wildcard_constraints:
        hap = '(1|2)'
    shell:
        'cat {input.hap} {input.un} > {output}'


rule canonical_dga_assemble_haplotypes_layout:
    input:
        fastq = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.{hap}.fastq.gz',
    output:
        layout = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/layout/{readset}.{hap}.ctg.lay.gz',
    wildcard_constraints:
        haplotype = '[h12\-un]+'
    log:
        'log/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.{hap}.layout.log',
    benchmark:
        'run/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/{readset}.{hap}.layout.rsrc',
    threads: 64
    params:
        param_preset = lambda wildcards: config['wtdbg2_presets'][wildcards.readset.rsplit('_', 1)[0]],
        out_prefix = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/layout/{readset}.{hap}'
    shell:
        'wtdbg2 -x {params.param_preset} -i {input.fastq} -g3g -t {threads} -o {params.out_prefix} &> {log}'


rule canonical_dga_assemble_haplotypes_consensus:
    input:
        layout = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/layout/{readset}.{hap}.ctg.lay.gz',
    output:
        'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/consensus/{readset}.{hap}.fasta',
    wildcard_constraints:
        hap = '[h12\-un]+'
    log:
        'log/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/consensus/{readset}.{hap}.log',
    benchmark:
        'run/output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{readset}/consensus/{readset}.{hap}.rsrc',
    threads: 64
    shell:
        'wtpoa-cns -t {threads} -i {input.layout} -o {output} &> {log}'