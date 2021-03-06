
def tech_comparison_determine_targets(wildcards):
    
    import sys

    module_outputs = {
        'assembly_delta': 'output/evaluation/HiFi_vs_CLR/{sample}_HiFi-{hap1}_vs_CLR-{hap2}.delta.gz',
        'assembly_diff': 'output/evaluation/HiFi_vs_CLR/{sample}_HiFi-{hap1}_vs_CLR-{hap2}.diff',
        'trio_delta_clr': 'output/evaluation/Trio_vs_CLR/{sample}_Trio-{hap1}_vs_CLR-{hap2}.delta.gz',
        'trio_delta_hifi': 'output/evaluation/Trio_vs_HiFi/{sample}_Trio-{hap1}_vs_HiFi-{hap2}.delta.gz',
        'trio_diff_clr': 'output/evaluation/Trio_vs_CLR/{sample}_Trio-{hap1}_vs_CLR-{hap2}.diff',
        'trio_diff_hifi': 'output/evaluation/Trio_vs_HiFi/{sample}_Trio-{hap1}_vs_HiFi-{hap2}.diff',
    }

    fix_wildcards = {
    }

    compute_results = set()

    # why?
    # this will fail if haplotype is unrecognized
    other_hap = {
        'h1-un': 'h2-un',
        'h2-un': 'h1-un'
    }

    search_path = os.path.join(os.getcwd(), 'output/evaluation/phased_assemblies')
    if not os.path.isdir(search_path):
        sys.stderr.write('\nNo phased assemblies at: {}\n'.format(search_path))
        return []
    for ps_assm in os.listdir(search_path):
        if ps_assm.startswith('v1'):
            version, new_name = ps_assm.split('_', 1)
            os.rename(os.path.join(search_path, ps_assm), os.path.join(search_path, new_name))
            assm_file = new_name
        else:
            assm_file = ps_assm
        if not assm_file.endswith('.fasta'):
            continue
        assm_base, hap, polisher, ext = assm_file.split('.')
        sample, assm_reads = assm_base.split('_', 1)

        if sample not in ['HG00733', 'NA19240', 'HG00514']:
            continue

        tmp = dict(fix_wildcards)
        tmp['sample'] = sample
        tmp['hap1'] = hap
        tmp['hap2'] = hap
        for target in module_outputs.values():
            if 'Trio' in target:
                tmp['hap1'] = 'hapA'
                fmt_target = target.format(**tmp)
                compute_results.add(fmt_target)
                tmp['hap1'] = 'hapB'
                fmt_target = target.format(**tmp)
                compute_results.add(fmt_target)
                tmp['hap1'] = hap
            else:
                fmt_target = target.format(**tmp)
                compute_results.add(fmt_target)

        tmp['hap2'] = other_hap[hap]
        for target in module_outputs.values():
            if 'Trio' in target:
                tmp['hap1'] = 'hapA'
                fmt_target = target.format(**tmp)
                compute_results.add(fmt_target)
                tmp['hap1'] = 'hapB'
                fmt_target = target.format(**tmp)
                compute_results.add(fmt_target)
                tmp['hap1'] = hap
            else:
                fmt_target = target.format(**tmp)
                compute_results.add(fmt_target)

    return sorted(compute_results)


localrules: master_tech_comparison


rule master_tech_comparison:
    input:
        tech_comparison_determine_targets


rule compute_assembly_delta:
    input:
        assm_ref = 'output/evaluation/phased_assemblies/{sample}_hgsvc_pbsq2-ccs_1000-pereg.{hap1}.racon-p2.fasta',
        assm_query = 'output/evaluation/phased_assemblies/{sample}_hgsvc_pbsq2-clr_1000-flye.{hap2}.arrow-p1.fasta',
    output:
        'output/evaluation/HiFi_vs_CLR/{sample}_HiFi-{hap1}_vs_CLR-{hap2}.delta'
    log:
        'log/output/evaluation/HiFi_vs_CLR/{sample}_HiFi-{hap1}_vs_CLR-{hap2}.log'
    benchmark:
        'run/output/evaluation/HiFi_vs_CLR/{sample}_HiFi-{hap1}_vs_CLR-{hap2}' + '.t{}.rsrc'.format(config['num_cpu_high'])
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads:
        config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 49152 + 24576 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: int((49152 + 24576 * attempt) / config['num_cpu_high']),
        runtime_hrs = lambda wildcards, attempt: 3 * attempt
    shell:
        'nucmer --maxmatch -l 100 -c 500 --threads={threads} --delta={output} '
            ' {input.assm_ref} {input.assm_query} &> {log}'


rule compute_trio_assembly_delta_CLR:
    input:
        assm_ref = 'references/assemblies/{sample}_hgsvc_pbil-trio_{hap1}.fasta',
        assm_query = 'output/evaluation/phased_assemblies/{sample}_hgsvc_pbsq2-clr_1000-flye.{hap2}.arrow-p1.fasta',
    output:
        'output/evaluation/Trio_vs_CLR/{sample}_Trio-{hap1}_vs_CLR-{hap2}.delta'
    log:
        'log/output/evaluation/Trio_vs_CLR/{sample}_Trio-{hap1}_vs_CLR-{hap2}.log'
    benchmark:
        'run/output/evaluation/Trio_vs_CLR/{sample}_Trio-{hap1}_vs_CLR-{hap2}' + '.t{}.rsrc'.format(config['num_cpu_high'])
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads:
        config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 262144 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: int(262144 * attempt / config['num_cpu_high']),
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    shell:
        'nucmer --maxmatch -l 100 -c 500 --threads={threads} --delta={output} '
            ' {input.assm_ref} {input.assm_query} &> {log}'


rule compute_trio_assembly_delta_HiFi:
    input:
        assm_ref = 'references/assemblies/{sample}_hgsvc_pbil-trio_{hap1}.fasta',
        assm_query = 'output/evaluation/phased_assemblies/{sample}_hgsvc_pbsq2-ccs_1000-pereg.{hap2}.racon-p2.fasta',
    output:
        'output/evaluation/Trio_vs_HiFi/{sample}_Trio-{hap1}_vs_HiFi-{hap2}.delta'
    log:
        'log/output/evaluation/Trio_vs_HiFi/{sample}_Trio-{hap1}_vs_HiFi-{hap2}.log'
    benchmark:
        'run/output/evaluation/Trio_vs_HiFi/{sample}_Trio-{hap1}_vs_HiFi-{hap2}' + '.t{}.rsrc'.format(config['num_cpu_high'])
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads:
        config['num_cpu_high']
    resources:
        mem_total_mb = lambda wildcards, attempt: 49152 + 24576 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: int((49152 + 24576 * attempt) / config['num_cpu_high']),
        runtime_hrs = lambda wildcards, attempt: 6 * attempt
    shell:
        'nucmer --maxmatch -l 100 -c 500 --threads={threads} --delta={output} '
            ' {input.assm_ref} {input.assm_query} &> {log}'


rule run_delta_diff:
    input:
        'output/evaluation/{tech1}_vs_{tech2}/{sample}_{tech1}-{hap1}_vs_{tech2}-{hap2}.delta'
    output:
        touch('output/evaluation/{tech1}_vs_{tech2}/{sample}_{tech1}-{hap1}_vs_{tech2}-{hap2}.diff')
    log:
        'log/output/evaluation/{tech1}_vs_{tech2}/{sample}_{tech1}-{hap1}_vs_{tech2}-{hap2}.diff.log'
    benchmark:
        'run/output/evaluation/{tech1}_vs_{tech2}/{sample}_{tech1}-{hap1}_vs_{tech2}-{hap2}.diff.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    priority: 1000
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt if attempt < 2 else attempt * 80
    params:
        out_dir = lambda wildcards, output: output[0].rsplit('.', 1)[0],
        out_prefix = lambda wildcards, output: os.path.join(
            output[0].rsplit('.', 1)[0],
            '{}_{}-{}_vs_{}-{}'.format(
                wildcards.sample,
                wildcards.tech1,
                wildcards.hap1,
                wildcards.tech2,
                wildcards.hap2
                )
            )
    shell:
        'mkdir -p {params.out_dir} && dnadiff -d {input} -p {params.out_prefix} &> {log}'