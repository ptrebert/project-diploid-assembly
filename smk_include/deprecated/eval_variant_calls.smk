
localrules: master_eval_variant_calls

rule master_eval_variant_calls:
    input:


rule whatshap_compare_variant_phasings_by_caller:
    """
    hap_reads = FASTQ file to be used for haplotype reconstruction
    dataset_info = depending on approach, combination of reference, vc_reads, and sseq_reads
    """
    input:
        phasing1 = 'output/diploid_assembly/{approach}/{var_caller1}_GQ{gq}_DP{dp}/{dataset_info}/{hap_reads}.phased.vcf',
        phasing2 = 'output/diploid_assembly/{approach}/{var_caller2}_GQ{gq}_DP{dp}/{dataset_info}/{hap_reads}.phased.vcf'
    output:
        comp = 'output/evaluation/variant_phasings/by_caller/{approach}/{var_caller1}_vs_{var_caller2}/setting_GQ{gq}_DP{dp}/{dataset_info}/{hap_reads}_{var_caller1}_vs_{var_caller2}.comp.tsv',
        swerr = 'output/evaluation/variant_phasings/by_caller/{approach}/{var_caller1}_vs_{var_caller2}/setting_GQ{gq}_DP{dp}/{dataset_info}/{hap_reads}_{var_caller1}_vs_{var_caller2}.swerr.bed',
        block = 'output/evaluation/variant_phasings/by_caller/{approach}/{var_caller1}_vs_{var_caller2}/setting_GQ{gq}_DP{dp}/{dataset_info}/{hap_reads}_{var_caller1}_vs_{var_caller2}.lblock.tsv',
    log:
        'log/output/evaluation/variant_phasings/by_caller/{approach}/{var_caller1}_vs_{var_caller2}/setting_GQ{gq}_DP{dp}/{dataset_info}/{hap_reads}_{var_caller1}_vs_{var_caller2}.comp.log'
    benchmark:
        'run/output/evaluation/variant_phasings/by_caller/{approach}/{var_caller1}_vs_{var_caller2}/setting_GQ{gq}_DP{dp}/{dataset_info}/{hap_reads}_{var_caller1}_vs_{var_caller2}.comp.rsrc'
    priority: 200
    shell:
        'whatshap --debug compare --names {wildcards.var_caller1},{wildcards.var_caller2} ' \
        ' --tsv-pairwise {output.comp} --switch-error-bed {output.swerr} ' \
        ' --longest-block-tsv {output.block} {input.phasing1} {input.phasing2} &> {log}'


rule whatshap_compare_variant_phasings_by_approach:
    """
    vc_reads = FASTQ file used for variant calling relative to reference
    hap_reads = FASTQ file to be used for haplotype reconstruction
    sseq_reads = FASTQ file used for strand-seq phasing
    """
    input:
        canonical = 'output/diploid_assembly/canonical/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}.phased.vcf',
        strandseq = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sseq_reads}/{hap_reads}.phased.vcf'
    output:
        comp = 'output/evaluation/variant_phasings/by_approach/canonical_vs_strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}_canonical_vs_{sseq_reads}.comp.tsv',
        swerr = 'output/evaluation/variant_phasings/by_approach/canonical_vs_strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}_canonical_vs_{sseq_reads}.swerr.bed',
        block = 'output/evaluation/variant_phasings/by_approach/canonical_vs_strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}_canonical_vs_{sseq_reads}.lblock.tsv',
    log:
        'log/output/evaluation/variant_phasings/by_approach/canonical_vs_strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}_canonical_vs_{sseq_reads}.comp.log'
    benchmark:
        'run/output/evaluation/variant_phasings/by_approach/canonical_vs_strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{hap_reads}_canonical_vs_{sseq_reads}.comp.rsrc'
    priority: 200
    shell:
        'whatshap --debug compare --names {wildcards.var_caller}_canonical,{wildcards.var_caller}_strandseq ' \
        ' --tsv-pairwise {output.comp} --switch-error-bed {output.swerr} ' \
        ' --longest-block-tsv {output.block} {input.canonical} {input.strandseq} &> {log}'