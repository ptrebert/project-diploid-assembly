
rule meryl_count_reference_kmers:
    input:
        sequence = ancient('/gpfs/project/projects/medbioinf/data/references/{reference}.fasta')
    output:
        kmer_db = directory('/gpfs/project/projects/medbioinf/data/references/{reference}.k15.db'),
        rep_kmer = '/gpfs/project/projects/medbioinf/data/references/{reference}.k15.rep-grt09998.txt'
    benchmark:
        'rsrc/output/kmer_db/{reference}.meryl-ref-kmer.rsrc'
    conda:
        '../../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 12288 + 4096 * attempt,
        mem_total_gb = lambda wildcards, attempt: int((12288 + 4096 * attempt)/1024) - 2,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        "meryl count k=15 threads={threads} memory={resources.mem_total_gb} output {output.kmer_db} {input.sequence} "
            " && "
            "meryl print greater-than distinct=0.9998 {output.kmer_db} > {output.rep_kmer}"