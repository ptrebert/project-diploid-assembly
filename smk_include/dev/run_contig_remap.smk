
REMAP_CONFIG = {
    'ref_assembly': 'GRCh38_HGSVC2_noalt',
    'min_mapq': 60,
    'annotations': [
        '20200723_GRCh38_p13_regions.bed',
        '20200723_GRCh38_p13_unresolved-issues.bed',
        'GRCh38_cytobands.bed'
    ]
}


def contig_remap_determine_targets(wildcards):

    contig_remap_targets = {
        'hap_assm': 'output/statistics/contigs_to_ref_aln/evaluation/phased_assemblies/{assembly}_map-to_{known_ref}.mapq{mapq}.stats',
        'nhr_assm': 'output/statistics/contigs_to_ref_aln/evaluation/nhrclust_assemblies/{assembly}_map-to_{known_ref}.mapq{mapq}.stats',
        'region_intersect': 'output/evaluation/completeness/{assembly_type}/{annotation}_OVL_{assembly}.tsv'
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
        assembly_type = os.path.split(path)[-1]
        tmp = dict(fix_wildcards)
        tmp['assembly_type'] = assembly_type
        for assm in os.listdir(path):
            if assm.startswith('v1'):
                version, new_name = assm.split('_', 1)
                os.rename(os.path.join(path, assm), os.path.join(path, new_name))
                assm_file = new_name
            else:
                assm_file = assm
            if not assm_file.endswith('.fasta'):
                continue
            tmp['assembly'] = assm_file.rsplit('.', 1)[0]
            fmt_target = remap_trg.format(**tmp)
            compute_results.add(fmt_target)
            for annotation in REMAP_CONFIG['annotations']:
                tmp['annotation'] = annotation
                fmt_target = contig_remap_targets['region_intersect'].format(**tmp)
                compute_results.add(fmt_target)
    
    return sorted(compute_results)


rule filter_merge_contig_alignments:
    input:
        'output/alignments/contigs_to_reference/evaluation/{assembly_type}/{bedfile}.bed'
    output:
        'output/alignments/contigs_to_reference/evaluation/{assembly_type}/{bedfile}.q60mrg.bed'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'egrep "\s60\s" {input} | bedtools merge -i - > {output}'


rule intersect_contig_alignments_annotation:
    input:
        'references/annotation/{annotation}.4c.bed',
        'output/alignments/contigs_to_reference/evaluation/{assembly_type}/{contigs}.q60mrg.bed'
    output:
        'output/evaluation/completeness/{assembly_type}/{annotation}_OVL_{contigs}.tsv'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools intersect -wao -a {input[0]} -b {input[1]} > {output}'


rule master:
    input:
        contig_remap_determine_targets