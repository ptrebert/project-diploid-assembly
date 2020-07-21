
include: 'prep_custom_references.smk'

def determine_possible_computations(wildcards):
    
    import sys

    module_outputs = {
        'assembly_delta': 'output/evaluation/HiFi_vs_CLR/{sample}_HiFi-{hap1}_vs_CLR-{hap2}.delta'
    }

    fix_wildcards = {
    }

    compute_results = set()

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

        if sample not in ['HG00733', 'NA19240']:
            continue

        tmp = dict(fix_wildcards)
        tmp['sample'] = sample
        tmp['hap1'] = hap
        if 'h1' in hap:
            tmp['hap2'] = 'h2-un'
        elif 'h2' in hap:
            tmp['hap1']
        else:
            raise ValueError('Unrecognized haplotype: {}'.format(assm_file))
    
        for target in module_outputs.values():
            fmt_target = target.format(**tmp)
            compute_results.add(fmt_target)
    print(compute_results)
    return sorted(compute_results)


localrules: master_kmer_analysis


rule master_tech_comparison:
    input:
        determine_possible_computations


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
        mem_total_mb = lambda wildcards, attempt: 65536 + 32768 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: int((65536 + 32768 * attempt) / config['num_cpu_high']),
        runtime_hrs = lambda wildcards, attempt: 23 * attempt
    shell:
        'nucmer --maxmatch -l 100 -c 500 --threads={threads} --delta={output} '
            ' {input.assm_ref} {input.assm_query} &> {log}'
