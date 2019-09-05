
localrules: ex_nihilo_read_splitting, ex_nihilo_polishing, ex_nihilo_consensus

rule ex_nihilo_read_splitting:
    output:
        protected(
            expand('output/haplotype_partitioning/splitting/HG00514.sequel2.pb-ccs.tag-{haplotype}.fastq.gz',
                    haplotype=['h1', 'h2', 'un'])
        ),
        protected(
            expand('output/haplotype_partitioning/splitting/HG00514.sequel2.pb-clr.tag-{haplotype}-un.fastq.gz',
                    haplotype=['h1', 'h2'])
        ),
        protected(
            expand('output/haplotype_partitioning/splitting/HG00514.sequel2.pb-clr-25.tag-{haplotype}-un.fastq.gz',
                    haplotype=['h1', 'h2'])
        )

    run:
        print('Ex nihilo splitting')


rule ex_nihilo_polishing:
    output:
        protected(
            expand('output/assembly_polishing/arrow_polished_contigs/{contig_sample}.{haplotype}.arrow-pass{round}.ctg.fa',
                    contig_sample=['HG00514.sequel2.pb-clr-25_to_HG00514.sequel2.pb-clr-25'],
                    haplotype=['tag-h1-un'],
                    round=[1])
        ),
        protected(
            expand('output/assembly_polishing/read_contig_alignments/{contig_sample}.{haplotype}.arrow-pass{round}.bam',
                    contig_sample=['HG00514.sequel2.pb-clr-25_to_HG00514.sequel2.pb-clr-25'],
                    haplotype=['tag-h1-un', 'tag-h2-un'],
                    round=[1])
        ),
        protected(
            expand('output/assembly_polishing/read_contig_alignments/{contig_sample}.{haplotype}.arrow-pass{round}.bam.{index}',
                    contig_sample=['HG00514.sequel2.pb-clr-25_to_HG00514.sequel2.pb-clr-25'],
                    haplotype=['tag-h1-un', 'tag-h2-un'],
                    index=['pbi'],
                    round=[1])
        )
    run:
        print('Ex nihilo polishing')


rule ex_nihilo_consensus:
    output:
        protected(
            expand('output/haplotype_assembly/consensus/HG00514.sequel2.pb-ccs.tag-{haplotype}.ctg.fa',
                    haplotype=['h1', 'h2', 'h1-un', 'h2-un'])
        ),
        protected(
            expand('output/haplotype_assembly/consensus/HG00514.sequel2.pb-clr-25.tag-{haplotype}.ctg.fa',
                haplotype=['h1-un', 'h2-un'])
        )
    run:
        print('Ex nihilo consensus')
