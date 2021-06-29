
rule map_racon_contigs_to_reference:
    input:
        contigs = 'output/assembly_polishing/racon_polished_contigs/{sample}.ctg.fa',
        reference = '/MMCI/TM/scratch/marschal/HG00514-assembly/ref/GRCh38.p13.renamed.fa',
        reference_index = '/MMCI/TM/scratch/marschal/HG00514-assembly/ref/GRCh38.p13.renamed.fa.fai'
    output:
        temp('output/assembly_analysis/contig_to_reference_alignment/{sample}.bam')
    log:
        minimap = 'log/assembly_analysis/contig_to_reference_alignment/{sample}.minimap.log',
        samtools = 'log/assembly_analysis/contig_to_reference_alignment/{sample}.samtools.log'
    benchmark: 'rsrc/assembly_analysis/contig_to_reference_alignment/{sample}.rsrc'
    threads: 48
    params:
        param_preset = 'asm10'
    run:
        exec = 'minimap2 -a -t {threads}'
        exec += ' -x {params.param_preset}'
        exec += ' {input.reference} {input.contigs}'
        exec += ' 2> {log.minimap}'
        exec += ' | samtools view -S -b -o {output} - '
        exec += ' 2> {log.samtools}'
        shell(exec)


rule map_arrow_contigs_to_reference:
    input:
        contigs = 'output/assembly_polishing/arrow_polished_contigs/{sample}.ctg.fa',
        reference = '/MMCI/TM/scratch/marschal/HG00514-assembly/ref/GRCh38.p13.renamed.fa',
        reference_index = '/MMCI/TM/scratch/marschal/HG00514-assembly/ref/GRCh38.p13.renamed.fa.fai'
    output:
        temp('output/assembly_analysis/contig_to_reference_alignment/{sample}.bam')
    log:
        minimap = 'log/assembly_analysis/contig_to_reference_alignment/{sample}.minimap.log',
        samtools = 'log/assembly_analysis/contig_to_reference_alignment/{sample}.samtools.log'
    benchmark: 'rsrc/assembly_analysis/contig_to_reference_alignment/{sample}.rsrc'
    threads: 48
    params:
        param_preset = 'asm10'
    run:
        exec = 'minimap2 -a -t {threads}'
        exec += ' -x {params.param_preset}'
        exec += ' {input.reference} {input.contigs}'
        exec += ' 2> {log.minimap}'
        exec += ' | samtools view -S -b -o {output} - '
        exec += ' 2> {log.samtools}'
        shell(exec)


rule contig_phasing:
    input:
        contig_alignment = 'output/assembly_analysis/contig_to_reference_alignment/{sample}.sort.bam',
        contig_alignment_index = 'output/assembly_analysis/contig_to_reference_alignment/{sample}.sort.bam.bai',
        reference = '/MMCI/TM/scratch/marschal/HG00514-assembly/ref/GRCh38.p13.renamed.fa',
        reference_index = '/MMCI/TM/scratch/marschal/HG00514-assembly/ref/GRCh38.p13.renamed.fa.fai',
        variants = '/MMCI/TM/scratch/marschal/HG00514-assembly/final-het-svns/ccs/HG00514.vcf.gz',
        variant_index = '/MMCI/TM/scratch/marschal/HG00514-assembly/final-het-svns/ccs/HG00514.vcf.gz.tbi'
    output:
        'output/assembly_analysis/whatshap/contig_phasing/{sample}.vcf.gz'
    log: 'log/assembly_analysis/whatshap/contig_phasing/{sample}.log'
    benchmark: 'rsrc/assembly_analysis/whatshap/contig_phasing/{sample}.rsrc'
    resources:
        mem_mb = 50 * (1024 * 1024)
    run:
        exec = 'whatshap phase --ignore-read-groups'
        exec += ' --reference {input.reference}'
        exec += ' --output {output}'
        exec += ' {input.variants}'
        exec += ' {input.contig_alignment}'
        exec += ' &> {log}'
        shell(exec)


rule contig_phasing_comparison:
    input:
        contig_phasings = 'output/assembly_analysis/whatshap/contig_phasing/{sample}.vcf.gz',
        reference_phasings = '/MMCI/TM/scratch/marschal/HG00514-assembly/input/trio-phasing/SH032.all.vcf.gz',
        reference_phasings_index = '/MMCI/TM/scratch/marschal/HG00514-assembly/input/trio-phasing/SH032.all.vcf.gz.tbi'
    output:
        table = 'output/assembly_analysis/whatshap/phasing_comparison/{sample}.eval.tsv',
        switch_errors = 'output/assembly_analysis/whatshap/phasing_comparison/{sample}.switch_errors.bed'
    log: 'log/assembly_analysis/whatshap/phasing_comparison/{sample}.eval.log'
    benchmark: 'rsrc/assembly_analysis/whatshap/phasing_comparison/{sample}.rsrc'
    priority: 100
    run:
        exec = 'whatshap compare --sample HG00514 --names \"ASSEMBLY,TRIO\"'
        exec += ' --tsv-pairwise {output.table}'
        exec += ' --switch-error-bed {output.switch_errors}'
        exec += ' {input.contig_phasings}'
        exec += ' {input.reference_phasings}'
        exec += ' &> {log}'
        shell(exec)
