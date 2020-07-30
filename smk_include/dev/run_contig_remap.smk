
localrules: master_contig_remap

REMAP_CONFIG = {
    'ref_assembly': 'GRCh38_HGSVC2_noalt',
    'min_mapq': 60,
    'annotations': [
        '20200723_GRCh38_p13_regions',
        '20200723_GRCh38_p13_unresolved-issues',
        'GRCh38_cytobands'
    ]
}


def contig_remap_determine_targets(wildcards):

    contig_remap_targets = {
        'hap_assm': 'output/statistics/contigs_to_ref_aln/evaluation/phased_assemblies/{assembly}_map-to_{known_ref}.mapq{mapq}.stats',
        'nhr_assm': 'output/statistics/contigs_to_ref_aln/evaluation/nhrclust_assemblies/{assembly}_map-to_{known_ref}.mapq{mapq}.stats',
        'region_intersect': 'output/evaluation/completeness/{assembly_type}/{annotation}_OVL_{assembly}_map-to_{known_ref}.tsv'
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
    
    cov_targets = {
        'region_avg': 'output/evaluation/completeness/hap_read_coverage/{annotation}_AVG_{readset}_map-to_hg38_GCA_p13.{hap}.tab'
    }

    cov_path = 'output/evaluation/hap_read_coverage'

    for fname in os.listdir(cov_path):
        if fname.startswith('v1'):
            version, new_name = fname.split('_', 1)
            os.rename(os.path.join(cov_path, fname), os.path.join(cov_path, new_name))
            bigwig = new_name
        else:
            bigwig = fname
        mapping, hap, _ = bigwig.split('.')
        readset, know_ref = mapping.split('_map-to_')
        formatter = {
            'readset': readset,
            'hap': hap
        }
        for trg in cov_targets.values():
            for annotation in REMAP_CONFIG['annotations']:
                formatter['annotation'] = annotation
                fmt_target = trg.format(**formatter)
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
        'output/alignments/contigs_to_reference/evaluation/{assembly_type}/{contigs}_map-to_{known_ref}.q60mrg.bed'
    output:
        'output/evaluation/completeness/{assembly_type}/{annotation}_OVL_{contigs}_map-to_{known_ref}.tsv'
    conda:
        '../../environment/conda/conda_biotools.yml'
    shell:
        'bedtools intersect -wao -a {input[0]} -b {input[1]} > {output}'


rule haploid_read_coverage_annotation:
    input:
        'references/annotation/{annotation}.4c.bed',
        'output/evaluation/hap_read_coverage/{readset}_map-to_hg38_GCA_p13.{hap}.bigWig'
    output:
        'output/evaluation/completeness/hap_read_coverage/{annotation}_AVG_{readset}_map-to_hg38_GCA_p13.{hap}.tab'
    conda:
        '../../environment/conda/conda_evaltools.yml'
    resources:
        mem_total_mb = lambda wildcards, attempt: 8192 if attempt < 2 else 32768 * attempt,
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 if attempt < 2 else 32768 * attempt
    shell:
        'bigWigAverageOverBed {input[1]} {input[0]} {output}'


rule master_contig_remap:
    input:
        contig_remap_determine_targets