
include: 'run_alignments.smk'

rule polish_assembled_contigs_pass1:
    input:
        use_for_correct = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/{pol_reads}.{hap}.fastq.gz',
        alignments = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/alignments/{hap_reads}/{pol_reads}_map-to_{reference}.{hap}.racon-p1.psort.sam',
        target_seqs = 'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/consensus/{hap_reads}.{hap}.fasta'
    output:
        'output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/{hap_reads}.{hap}.racon-p1.fasta'
    log:
        'log/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/{hap_reads}.{hap}.racon-p1.log'
    benchmark:
        'run/output/diploid_assembly/strandseq/{var_caller}_GQ{gq}_DP{dp}/{reference}/{vc_reads}/{sts_reads}/polishing/{pol_reads}/{hap_reads}.{hap}.racon-p1.rsrc'
    wildcard_constraints:
        gq = '[0-9]+',
        dp = '[0-9]+',
        reference = '[\w\-]+',
        vc_reads = '[\w\-]+',
        sts_reads = '[\w\-]+',
        hap_reads = '[\w\-]+',
        pol_reads = '[\w\-]+',
        hap = '[h12un\-]+',
    threads: 12
    shell:
        'racon --threads {threads} --include-unpolished {input.use_for_correct} {input.alignments} {input.target_seqs} > {output} 2> {log}'

#
#rule polish_assembled_contigs_pass2:
#    input:
#        use_for_correct = 'output/haplotype_partitioning/splitting/{read_sample}.{haplotype}.fastq.gz',
#        alignments = 'output/assembly_polishing/read_contig_alignments/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass2.sam',
#        correct_target = 'output/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass1.ctg.fa'
#    output:
#        'output/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass2.ctg.fa'
#    log: 'log/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass2.log'
#    benchmark: 'run/assembly_polishing/racon_ploished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon-pass2.rsrc'
#    wildcard_constraints:
#        read_sample = '[A-Z0-9\.\-]+',
#        contig_sample = '[A-Z0-9\.\-]+'
#    threads: 12
#    run:
#        exec = ' /home/pebert/work/code/github/external/racon/build/bin/racon --threads {threads}'
#        exec += ' --include-unpolished'
#        exec += ' {input.use_for_correct}'
#        exec += ' {input.alignments}'
#        exec += ' {input.correct_target}'
#        exec += ' > {output}'
#        exec += ' 2> {log}'
#        shell(exec)
#
