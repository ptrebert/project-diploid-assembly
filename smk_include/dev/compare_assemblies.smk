
rule compute_delta_between_assemblies:
    input:
        contigs1 = 'output/haplotype_assembly/consensus/{assm1}.{haplotype}.ctg.fa',
        contigs2 = 'output/assembly_polishing/racon_polished_contigs/{assm2}.{haplotype}.racon-pass{round}.ctg.fa',
    output:
        delta = 'output/assembly_analysis/mummer/{assm1}_vs_{assm2}.racon-pass{round}.{haplotype}.delta'
    log: 'log/assembly_analysis/mummer/{assm1}_vs_{assm2}.racon-pass{round}.{haplotype}.log'
    benchmark: 'run/assembly_analysis/mummer/{assm1}_vs_{assm2}.racon-pass{round}.{haplotype}.rsrc'
    threads: 24
    run:
        exec = 'nucmer --maxmatch -l 100 -c 500'
        exec += ' --threads={threads}'
        exec += ' --delta={output}'
        exec += ' {input.contigs1} {input.contigs2}'
        exec += ' &> {log}'
        shell(exec)