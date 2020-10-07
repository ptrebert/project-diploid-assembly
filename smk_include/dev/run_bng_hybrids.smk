
localrules: master_bng_hybrids

HYBRID_CONFIG = {
    'chroms': ['chr' + str(i) for i in range(1, 23)] + ['chrXY', 'chrUn']
}


def bng_hybrids_determine_targets(wildcards):

    hybrid_targets = {
        'contig_stats': 'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid.contig-stats.tsv'
    }

    search_paths = {
        (os.path.join(os.getcwd(), 'output/evaluation/scaffolded_assemblies'), 'contig_stats'),
    }

    compute_results = set()

    for path, trg_type in search_paths:
        if not os.path.isdir(path):
            continue
        hybrid_trg = hybrid_targets[trg_type]
        for hybrid_file in os.listdir(path):
            if not hybrid_file.endswith('.bng-hybrid.agp'):
                continue
            tmp = dict()
            tmp['assembly'] = hybrid_file.rsplit('.', 2)[0]
            fmt_target = hybrid_trg.format(**tmp)
            compute_results.add(fmt_target)
            
    return sorted(compute_results)


rule summarize_hybrid_statistics:
    input:
        agp = 'output/evaluation/scaffolded_assemblies/{assembly}.bng-hybrid.agp',
        fasta = 'output/evaluation/scaffolded_assemblies/{assembly}.bng-scaffolds.fasta',
        discard = 'output/evaluation/scaffolded_assemblies/{assembly}.bng-unsupported.fasta',
        ctg_ref_aln = 'output/alignments/contigs_to_reference/evaluation/phased_assemblies/{assembly}_map-to_GRCh38_HGSVC2_noalt.bed'
    output:
        contig_stats = 'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid.contig-stats.tsv',
        layout = 'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid.scaffold-layout.tsv',
        scaffolds = 'output/evaluation/bng_hybrids/{assembly}/{assembly}.hybrid.scaffolds.wg.fasta',
        contig_fasta = expand(
            'output/evaluation/bng_hybrids/{{assembly}}/{{assembly}}.hybrid.contigs.{chrom}.fasta',
            chrom=HYBRID_CONFIG['chroms']
        )
    log:
        'log/output/evaluation/bng_hybrids/{assembly}.hybrid.stats.log',
    conda:
        '../../environment/conda/conda_evaltools.yml'
    resources:
        mem_per_cpu_mb = lambda wildcards, attempt: 8192 * attempt,
        mem_total_mb = lambda wildcards, attempt: 8192 * attempt,
    params:
        exec = lambda wildcards: find_script_path('process_bng_hybrid.py'),
        out_prefix = lambda wildcards, output: output.contig_stats.rsplit('.', 3)[0] 
    shell:
        '{params.exec} --agp-file {input.agp} --fasta-file {input.fasta} '
            '--bed-file {input.ctg_ref_aln} --output {params.out_prefix} &> {log}'


rule master_bng_hybrids:
    input:
        bng_hybrids_determine_targets
