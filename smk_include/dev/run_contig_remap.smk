
REMAP_CONFIG = {
    'ref_assembly': 'GRCh38_HGSVC2_noalt',
    'min_mapq': 60
}


def contig_remap_determine_targets(wildcards):

    contig_remap_targets = {
        'hap_assm': 'output/statistics/contigs_to_ref_aln/evaluation/phased_assemblies/{assembly}_map-to_{known_ref}.mapq{mapq}.stats',
        'nhr_assm': 'output/statistics/contigs_to_ref_aln/evaluation/nhrclust_assemblies/{assembly}_map-to_{known_ref}.mapq{mapq}.stats'
    }

    fix_wildcards = {
        'known_ref': REMAP_CONFIG['ref_assembly'],
        'mapq': REMAP_CONFIG['min_mapq']
    }

    compute_results = set()

    search_paths = [
        (os.path.join(os.getcwd(), 'output/evaluation/phased_assemblies'), 'hap_assm'),
        (os.path.join(os.getcwd(), 'output/evaluation/nhrclust_assemblies'), 'nhr_assm'),
    ]

    # limit to phased assemblies for now
    for path, trg_type in search_paths:
        remap_trg = contig_remap_targets[trg_type]
        tmp = dict(fix_wildcards)
        for assm in os.listdir(path):
            if assm.startswith('v1'):
                version, new_name = assm.split('_', 1)
                os.rename(os.path.join(search_path, assm), os.path.join(search_path, new_name))
                assm_file = new_name
            else:
                assm_file = assm
            if not assm_file.endswith('.fasta'):
                continue
            tmp['assembly'] = assm_file.rsplit('.', 1)[0]
            fmt_target = remap_trg.format(**tmp)
            compute_results.add(fmt_target)
    
    return sorted(compute_results)


rule master:
    input:
        contig_remap_determine_targets