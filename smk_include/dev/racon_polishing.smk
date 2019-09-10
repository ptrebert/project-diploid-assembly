
# For potential future use

rule polish_assembled_contigs_pass1:
    input:
        use_for_correct = 'output/haplotype_partitioning/splitting/{read_sample}.{haplotype}.fastq.gz',
        alignments = 'output/assembly_polishing/read_contig_alignments/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass1.sam',
        correct_target = 'output/haplotype_assembly/consensus/{contig_sample}.{haplotype}.ctg.fa'
    output:
        'output/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass1.ctg.fa'
    log: 'log/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass1.log'
    benchmark: 'run/assembly_polishing/racon_ploished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass1.rsrc'
    wildcard_constraints:
        read_sample = '[A-Z0-9\.\-]+'
    threads: 12
    run:
        exec = ' /home/pebert/work/code/github/external/racon/build/bin/racon --threads {threads}'
        exec += ' --include-unpolished'
        exec += ' {input.use_for_correct}'
        exec += ' {input.alignments}'
        exec += ' {input.correct_target}'
        exec += ' > {output}'
        exec += ' 2> {log}'
        shell(exec)


rule polish_assembled_contigs_pass2:
    input:
        use_for_correct = 'output/haplotype_partitioning/splitting/{read_sample}.{haplotype}.fastq.gz',
        alignments = 'output/assembly_polishing/read_contig_alignments/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass2.sam',
        correct_target = 'output/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass1.ctg.fa'
    output:
        'output/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass2.ctg.fa'
    log: 'log/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass2.log'
    benchmark: 'run/assembly_polishing/racon_ploished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass2.rsrc'
    wildcard_constraints:
        read_sample = '[A-Z0-9\.\-]+',
        contig_sample = '[A-Z0-9\.\-]+'
    threads: 12
    run:
        exec = ' /home/pebert/work/code/github/external/racon/build/bin/racon --threads {threads}'
        exec += ' --include-unpolished'
        exec += ' {input.use_for_correct}'
        exec += ' {input.alignments}'
        exec += ' {input.correct_target}'
        exec += ' > {output}'
        exec += ' 2> {log}'
        shell(exec)


rule polish_arrow_contigs_racon_pass1:
    input:
        use_for_correct = 'output/haplotype_partitioning/splitting/{read_sample}.{haplotype}.fastq.gz',
        alignments = 'output/assembly_polishing/read_contig_alignments/{read_sample}_to_{contig_sample}.{haplotype}.arrow-pass{round}.racon-pass{round}.sam',
        correct_target = 'output/assembly_polishing/arrow_polished_contigs/{contig_sample}.{haplotype}.arrow-pass{round}.ctg.fa'
    output:
        'output/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.arrow-pass{round}.racon-pass{round}.ctg.fa'
    log: 'log/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.arrow-pass{round}.racon-pass{round}.log'
    benchmark: 'run/assembly_polishing/racon_ploished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.arrow-pass{round}.racon-pass{round}.rsrc'
    wildcard_constraints:
        read_sample = '[A-Za-z0-9\.\-]+',
    threads: 12
    run:
        exec = ' /home/pebert/work/code/github/external/racon/build/bin/racon --threads {threads}'
        exec += ' --include-unpolished'
        exec += ' {input.use_for_correct}'
        exec += ' {input.alignments}'
        exec += ' {input.correct_target}'
        exec += ' > {output}'
        exec += ' 2> {log}'
        shell(exec)
