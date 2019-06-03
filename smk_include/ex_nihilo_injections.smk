
rule ex_nihilo_read_splitting:
    output:
        expand('output/haplotype_partitioning/splitting/HG00514.sequel2.pb-ccs.tag-{haplotype}.fastq.gz',
                haplotype=['h1', 'h2', 'un']),
        expand('output/haplotype_partitioning/splitting/HG00514.sequel2.pb-clr.tag-{haplotype}-un.fastq.gz',
                haplotype=['h1', 'h2'])
    run:
        print('Ex nihilo splitting')


rule ex_nihilo_consensus:
    output:
        expand('output/haplotype_assembly/consensus/HG00514.sequel2.pb-ccs.tag-{haplotype}.ctg.fa',
                haplotype=['h1', 'h2', 'h1-un', 'h2-un'])
    run:
        print('Ex nihilo consensus')
