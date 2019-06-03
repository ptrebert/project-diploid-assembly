
import collections as col

REFKEY = config['use_ref_genome']
REFNAME = REFKEY + '.fa'
REFCHROM = config['chromosomes'][REFKEY]

SAMPLES = config['samples']

# this needs to be sorted out when finalizing the pipeline
ruleorder: merge_phased_variants > prepare_trio_phased_variants
ruleorder: merge_chunk_part_read_files > merge_partitioned_read_files

include: 'smk_include/aux_utilities.smk'
include: 'smk_include/ex_nihilo_injections.smk'
include: 'smk_include/dl_reference_data.smk'
include: 'smk_include/dl_read_data_chunks.smk'
include: 'smk_include/dl_read_data_parts.smk'

rule master:
    input:
        # input validation
        expand('output/fastq_validation/{sample}.stats.pck',
                sample=SAMPLES),

        # assembled haplotypes
        # - only haplotypes
        # - haplotypes plus untagged reads
#        expand('output/haplotype_assembly/consensus/{sample}.tag-{haplotype}{untagged}.ctg.fa',
#               sample=SAMPLES,
#               haplotype=['h1', 'h2'],
#               untagged=['', '-un']),

#        expand('output/haplotype_assembly/consensus/{sample}.tag-{haplotype}{untagged}.ctg.fa',
#               sample=['HG00514.sequel2.pb-clr'],
#               haplotype=['h1', 'h2'],
#               untagged=['-un']),
        expand('output/assembly_analysis/quastlg/{sample}.tag-{haplotype}{untagged}/report.pdf',
               sample=['HG00514.sequel2.pb-clr'],
               haplotype=['h1', 'h2'],
               untagged=['-un']),

        expand('output/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon.ctg.fa',
                read_sample=['HG00514.sequel2.pb-ccs'],
                contig_sample=['HG00514.sequel2.pb-clr'],
                haplotype=['tag-h1-un', 'tag-h2-un']),

        expand('output/assembly_analysis/quastlg/{read_sample}_to_{contig_sample}.{haplotype}.racon/report.pdf',
               contig_sample=['HG00514.sequel2.pb-clr'],
               read_sample=['HG00514.sequel2.pb-ccs'],
               haplotype=['tag-h1-un', 'tag-h2-un']),

        # assembled haplotypes
        # - only haplotypes
        # - haplotypes plus untagged reads
#        expand('output/assembly_analysis/mummer/{sample}.tag-{haplotype}{untagged}.coord.tab',
#               sample=SAMPLES,
#               haplotype=['h1', 'h2'],
#               untagged=['', '-un'])
#        expand('output/assembly_analysis/quastlg/{sample}.tag-{haplotype}{untagged}/report.pdf',
#               sample=SAMPLES,
#               haplotype=['h1', 'h2'],
#               untagged=['', '-un']),
#        expand('output/assembly_analysis/quastlg/{sample}.tag-{haplotype}{untagged}/report.pdf',
#                sample=['HG00514.sequel2.pb-ccs'],
#                haplotype=['h1', 'h2'],
#                untagged=['', '-un']),
#
#        expand('output/assembly_analysis/mummer/{sample}.tag-{haplotype}{untagged}.coord.tab',
#               sample=['HG00514.sequel2.pb-ccs'],
#               haplotype=['h1', 'h2'],
#               untagged=['', '-un']),
#        expand('output/fastq_validation/{sample}.tag-{haplotype}.stats.pck',
#                sample=['HG00514.sequel2.pb-ccs'],
#                haplotype=['h1', 'h2', 'un']),

#
#        # Phase variants
#        expand('output/phased_variants/{sample}.{project}.{chromosome}.vcf.gz',
#                sample=['HG00733'], project=['pangen'],
#                chromosome=REFCHROM),
#
#        # Merge phased variants
#        expand('output/merged_phased_variants/{sample}.{project}.vcf.gz',
#                sample=['HG00733'], project=['pangen']),
#        expand('output/merged_phased_variants/{sample}.{project}.vcf.gz',
#                sample=['HG002'], project=['pangen', 'ucsc1']),
#
#        # Split by haplotype
#        expand('output/haplosplit_reads/{sample}.{project}.ont-ul.tag-{tag}.fastq.gz',
#                sample=['HG00733'], project=['pangen'], tag=['h1', 'h2', 'un']),
#        # GIAB sample still buggy
#        expand('output/haplosplit_reads/{sample}.{project}.ont-ul.tag-{tag}.fastq.gz',
#                sample=['HG002'], project=['pangen', 'ucsc1'],
#                tag=['h1', 'h2', 'un']),
#
#        # Compute MD5 before uploading
#        expand('output/haplosplit_reads/{sample}.{project}.ont-ul.tag-{tag}.fastq.gz.md5',
#                sample=['HG00733'], project=['pangen'], tag=['h1', 'h2', 'un']),
#        expand('output/haplosplit_reads/{sample}.{project}.ont-ul.tag-{tag}.fastq.gz.md5',
#                sample=['HG002'], project=['pangen', 'ucsc1'],
#                tag=['h1', 'h2', 'un']),
#
#         expand('output/assembled_haplotypes/{sample}.{project}.ont-ul.tag-{haplotype}.ctg.lay.gz',
#                sample=['HG002'], project=['pangen'],
#                haplotype=['h1', 'h2']),
##         expand('output/assembled_haplotypes/{sample}.{project}.ont-ul.tag-{haplotype}.ctg.lay.gz',
##                sample=['HG00733'], project=['pangen'],
##                haplotype=['h1', 'h2']),
#
#         expand('output/assembled_haplotypes/{sample}.{project}.ont-ul.tag-{haplotype}-un.ctg.lay.gz',
#                sample=['HG002'], project=['pangen'],
#                haplotype=['h1', 'h2']),
##         expand('output/assembled_haplotypes/{sample}.{project}.ont-ul.tag-{haplotype}-un.ctg.lay.gz',
##                sample=['HG00733'], project=['pangen'],
##                haplotype=['h1', 'h2'])
#

#
#         expand('output/assembly_analysis/mummer/{sample}.{project}.ont-ul.tag-{haplotype}{untagged}.delta',
#                sample=['HG002'], project=['pangen'],
#                haplotype=['h1', 'h2'], untagged=['', '-un']),
#
#         expand('output/assembly_analysis/mummer/{sample}.{project}.ont-ul.tag-{haplotype}{untagged}.delta.gz',
#                sample=['HG002'], project=['pangen'],
#                haplotype=['h1', 'h2'], untagged=['', '-un']),
#
#         expand('output/assembly_analysis/mummer/{sample}.{project}.ont-ul.tag-{haplotype}{untagged}.coord.tab',
#                sample=['HG002'], project=['pangen'],
#                haplotype=['h1', 'h2'], untagged=['', '-un']),
#         expand('output/assembly_analysis/mummer/{sample}.{project}.ont-ul.tag-{haplotype}{untagged}.coord.tab',
#                sample=['HG00733'], project=['pangen'],
#                haplotype=['h1', 'h2'], untagged=['', '-un'])

    message: 'Executing ALL'


onsuccess:
    body_text = "Nothing to report\n{log}"
    if config['notify']:
        shell('mail -s "[Snakemake] diploid genome assembly - success" {} <<< "{}"'.format(config['notify_email'], body_text))

onerror:
    if config['notify']:
        shell('mail -s "[Snakemake] diploid genome assembly - ERRROR" {} < {{log}}'.format(config['notify_email']))


rule prepare_trio_phased_variants:
    input:
        vcf = 'references/trio-phased-variants/AJ.vcf.gz',
        tbi = 'references/trio-phased-variants/AJ.vcf.gz.tbi'
    output:
        vcf = 'output/merged_phased_variants/{sample}.{project}.vcf.gz',
        tbi = 'output/merged_phased_variants/{sample}.{project}.vcf.gz.tbi',
    log: 'log/prep_trio_phased/{sample}.{project}.log'
    run:
        exec = 'bcftools view -o {output.vcf} -O z --samples {wildcards.sample} {input.vcf} &> {log}'
        exec += ' && '
        exec += 'bcftools index --tbi --output-file {output.tbi} {output.vcf} &>> {log}'
        shell(exec)


####################################################
# The following rules represent the actual pipeline
# Everything above is subject to change depending
# on the input data and final phasing strategy
####################################################


rule map_reads_to_reference:
    input:
        reference = 'references/' + REFNAME,
        read_data = 'input/read_data/diploid_assembly_input/{sample}.fastq.gz'
    output:
        temp('output/temp/alignments/{sample}.cram')
    log:
        minimap = 'log/alignments/{sample}.minimap.log',
        samtools = 'log/alignments/{sample}.samtools.log'
    benchmark: 'run/alignments/{sample}.minimap.samtools.rsrc'
    params:
        cache_path = os.path.join(os.getcwd(), config['ref_cache_folder']),
        param_preset = lambda wildcards: config['minimap2_presets'][wildcards.sample]
    threads: 48
    run:
        exec = 'REF_CACHE={params.cache_path}'
        exec += ' minimap2 -t {threads}'
        exec += ' -R \'@RG\\tID:1\\tSM:{wildcards.sample}\''
        exec += ' -a'  # output in SAM format
        exec += ' -x {params.param_preset}'
        exec += ' {input.reference}'
        exec += ' {input.read_data} 2> {log.minimap}'
        exec += ' | samtools view'
        exec += ' -T {input.reference}'
        exec += ' -C -'
        exec += ' > {output}'
        exec += ' 2> {log.samtools}'
        shell(exec)


rule map_haplotype_reads_to_contigs:
    input:
        contigs = 'output/haplotype_assembly/consensus/{contig_sample}.{haplotype}.ctg.fa',
        read_data = 'output/haplotype_partitioning/splitting/{read_sample}.{haplotype}.fastq.gz'
    output:
        'output/assembly_polishing/read_contig_alignments/{read_sample}_to_{contig_sample}.{haplotype}.paf'
    log: 'log/assembly_polishing/read_contig_alignments/{read_sample}_to_{contig_sample}.{haplotype}.aln.log'
    benchmark: 'run/assembly_polishing/read_contig_alignments/{read_sample}_to_{contig_sample}.{haplotype}.rsrc'
    wildcard_constraints:
        hap = 'h[12]'
    params:
        param_preset = lambda wildcards: config['minimap2_presets'][wildcards.read_sample]
    threads: 48
    run:
        exec = ' minimap2 -t {threads}'
        exec += ' -x {params.param_preset}'
        exec += ' {input.contigs}'
        exec += ' {input.read_data}'
        exec += ' > {output} 2> {log}'
        shell(exec)


# sorting takes ~1 hour / 50 GiB with 20 cores and 5GiB / core
# ATTENTION: current settings result in 20 * 20 GiB RAM => 400 GiB required!
rule sort_cram_alignments:
    input:
        'output/temp/alignments/{sample}.cram'
    output:
        'output/alignments/{sample}.sort.cram'
    log: 'log/aln_sort/{sample}.sort.log'
    benchmark: 'run/alignments/{sample}.sort.rsrc'
    threads: 20
    resources:
        mem_mb = 20 * 20480
    params:
        mem_limit = '20G',
        cache_path = os.path.join(os.getcwd(), config['ref_cache_folder'])
    run:
        exec = 'REF_CACHE={params.cache_path} '
        exec += 'samtools sort -o {output} -m {params.mem_limit}'
        exec += ' -@ {threads} {input} &> {log}'
        shell(exec)


rule phase_variants:
    input:
        vcf = 'output/variant_calls/{sample}/{sample}.{chromosome}.vcf.gz',
        tbi = 'output/variant_calls/{sample}/{sample}.{chromosome}.vcf.gz.tbi',
        cram = 'output/alignments/{sample}.sort.cram',
        crai = 'output/alignments/{sample}.sort.cram.crai',
        ref = 'references/' + REFNAME,
    output:
        vcf = 'output/phased_variants/{sample}.{chromosome}.vcf.gz',
    log: 'log/phase_variants/{sample}.{chromosome}.vcf.log'
    run:
        exec = 'whatshap phase --reference {input.ref}'
        exec += ' --chromosome {wildcards.chromosome}'
        exec += ' --output {output.vcf}'
        exec += ' {input.vcf} {input.cram}'
        exec += ' &> {log}'
        shell(exec)


rule merge_phased_variants:
    input:
        expand('output/phased_variants/{{sample}}.{chromosome}.vcf.gz',
                chromosome=REFCHROM)
    output:
        vcf = 'output/merged_phased_variants/{sample}.{project}.vcf.gz',
        tbi = 'output/merged_phased_variants/{sample}.{project}.vcf.gz.tbi'
    log: 'log/merge_phased/{sample}.{project}.log'
    run:
        exec = 'bcftools concat -o {output.vcf} -O z {input} &> {log}'
        exec += ' && '
        exec += 'bcftools index --tbi --output-file {output.tbi} {output.vcf} &>> {log}'
        shell(exec)


rule haplotype_partitioning_read_tagging:
    input:
        vcf = 'output/merged_phased_variants/{sample}.vcf.gz',
        tbi = 'output/merged_phased_variants/{sample}.vcf.gz.tbi',
        cram = 'output/alignments/{sample}.sort.cram',
        crai = 'output/alignments/{sample}.sort.cram.crai',
        ref = 'references/' + REFNAME,
    output:
        cram = 'output/haplotype_partitioning/tagging/{sample}.haplotag.cram',
        taglist = 'output/haplotype_partitioning/tagging/{sample}.haplotag.tsv.gz'
    log: 'log/haplotag/{sample}.haplotag.log'
    conda:
        "environment/conda/wh_split.yml"
    shell:
        "whatshap --debug haplotag --output {output.cram} --reference {input.ref} --output-haplotag-list {output.taglist} {input.vcf} {input.cram} &> {log}"


rule haplotype_partitioning_read_splitting:
    input:
        fastq = 'input/read_data/diploid_assembly_input/{sample}.fastq.gz',
        taglist = 'output/haplotype_partitioning/tagging/{sample}.haplotag.tsv.gz'
    output:
        haplo1 = 'output/haplotype_partitioning/splitting/{sample}.tag-h1.fastq.gz',
        haplo2 = 'output/haplotype_partitioning/splitting/{sample}.tag-h2.fastq.gz',
        untag = 'output/haplotype_partitioning/splitting/{sample}.tag-un.fastq.gz'
    log:
        'log/haplosplit/{sample}.haplosplit.log'
    conda:
        "environment/conda/wh_split.yml"
    shell:
        "whatshap --debug split --pigz --output-h1 {output.haplo1} --output-h2 {output.haplo2} --output-untagged {output.untag} {input.fastq} {input.taglist} &> {log}"


rule haplotype_assembly_haplotypes_only:
    input:
        fastq = 'output/haplotype_partitioning/splitting/{sample}.tag-{haplotype}.fastq.gz',
    output:
        layout = 'output/haplotype_assembly/assembly/{sample}.tag-{haplotype}/{sample}.tag-{haplotype}.ctg.lay.gz',
    wildcard_constraints:
        sample = '[A-Za-z0-9_\-\.]+',
        haplotype = '[h12\-un]+'
    log: 'log/haplotype_assembly/{sample}.tag-{haplotype}.layout.log'
    benchmark: 'run/haplotype_assembly/layout/{sample}.tag-{haplotype}.layout.rsrc'
    threads: 48
    params:
        param_preset = lambda wildcards: config['wtdbg2_presets'][wildcards.sample]
    run:
        exec = 'wtdbg2 -x {params.param_preset}'
        exec += ' -i {input.fastq}'
        exec += ' -g3g -t {threads}'  # approx genome size
        exec += ' -o output/haplotype_assembly/assembly/{wildcards.sample}.tag-{wildcards.haplotype}/{wildcards.sample}.tag-{wildcards.haplotype}'
        exec += ' &> {log}'
        shell(exec)


#rule haplotype_assembly_haplotypes_plus_untagged:
#    input:
#        hap_fastq = 'output/haplotype_partitioning/splitting/{sample}.tag-{haplotype}.fastq.gz',
#        un_fastq = 'output/haplotype_partitioning/splitting/{sample}.tag-un.fastq.gz',
#    output:
#        layout = 'output/haplotype_assembly/assembly/{sample}.tag-{haplotype}-un/{sample}.tag-{haplotype}-un.ctg.lay.gz',
#    wildcard_constraints:
#        sample = '[A-Za-z0-9_\-\.]+',
#        haplotype = '[h12]+'
#    log: 'log/haplotype_assembly/layout/{sample}.tag-{haplotype}-un.layout.log'
#    benchmark: 'run/haplotype_assembly/layout/{sample}.tag-{haplotype}-un.layout.rsrc'
#    threads: 48
#    params:
#        param_preset = lambda wildcards: config['wtdbg2_presets'][wildcards.sample]
#    run:
#        exec = 'wtdbg2 -x {params.param_preset}'
#        exec += ' -i {input.hap_fastq} -i {input.un_fastq}'
#        exec += ' -g3g -t {threads}'  # approx genome size
#        exec += ' -o output/haplotype_assembly/assembly/{wildcards.sample}.tag-{wildcards.haplotype}-un/{wildcards.sample}.tag-{wildcards.haplotype}-un'
#        exec += ' &> {log}'
#        shell(exec)


rule derive_assembly_consensus:
    input:
        ctg_layout = 'output/haplotype_assembly/assembly/{hapassm}/{hapassm}.ctg.lay.gz'
    output:
        'output/haplotype_assembly/consensus/{hapassm}.ctg.fa'
    wildcard_constraints:
        hapassm = '[A-Za-z0-9\.\-]+'
    log: 'log/haplotype_consensus/{hapassm}.log'
    benchmark: 'run/haplotype_consensus/{hapassm}.rsrc'
    threads: 48
    run:
        exec = 'wtpoa-cns -t {threads}'
        exec += ' -i {input.ctg_layout}'
        exec += ' -o {output} &> {log}'
        shell(exec)


rule compute_assembly_delta:
    input:
        contigs = 'output/haplotype_assembly/consensus/{hapassm}.ctg.fa',
        reference = 'references/' + REFNAME
    output:
        delta = 'output/assembly_analysis/mummer/{hapassm}.delta'
    log: 'log/assembly_analysis/mummer/{hapassm}.log'
    benchmark: 'run/assembly_analysis/mummer/{hapassm}.rsrc'
    threads: 24
    run:
        exec = 'nucmer --maxmatch -l 100 -c 500'
        exec += ' --threads={threads}'
        exec += ' {input.reference} {input.contigs}'
        exec += ' --delta={output}'
        exec += ' &> {log}'
        shell(exec)


rule convert_delta_to_coordinates:
    input:
        delta = 'output/assembly_analysis/mummer/{hapassm}.delta'
    output:
        tab = 'output/assembly_analysis/mummer/{hapassm}.coord.tab'
    params:
        minlen = 10000
    run:
        exec = 'show-coords -l'  # incl. seq. len. info.
        exec += ' -L {params.minlen}'  # min. aln. len.
        exec += ' -r'  # sort by reference ID
        exec += ' -T'  # switch to tab-sep output
        exec += ' {input.delta} > {output.tab}'
        shell(exec)


rule quast_assembly_statistics:
    input:
        contigs = 'output/haplotype_assembly/consensus/{hapassm}.ctg.fa',
        reference = 'references/' + REFNAME,
        genes = 'references/gencode.v30.basic.annotation.gff3.gz'
    output:
        pdf_report = 'output/assembly_analysis/quastlg/{hapassm}/report.pdf',
        html_icarus = 'output/assembly_analysis/quastlg/{hapassm}/icarus.html'
    threads: 16
    log: 'log/assembly_analysis/quastlg/{hapassm}.log'
    benchmark: 'run/assembly_analysis/quastlg/{hapassm}.rsrc'
    priority: 100
    params:
        output_dir = 'output/assembly_analysis/quastlg/{hapassm}'
    run:
        exec = 'quast-lg.py --threads {threads}'
        exec += ' -r {input.reference}'
        exec += ' --features gene:{input.genes}'
        exec += ' --output-dir {params.output_dir}'
        exec += ' {input.contigs}'
        exec += ' &> {log}'
        shell(exec)


rule polish_assembled_contigs:
    input:
        use_for_correct = 'output/haplotype_partitioning/splitting/{read_sample}.{haplotype}.fastq.gz',
        alignments = 'output/assembly_polishing/read_contig_alignments/{read_sample}_to_{contig_sample}.{haplotype}.paf',
        correct_target = 'output/haplotype_assembly/consensus/{contig_sample}.{haplotype}.ctg.fa'
    output:
        'output/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon.ctg.fa'
    log: 'log/assembly_polishing/racon_polished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon.log'
    benchmark: 'run/assembly_polishing/racon_ploished_contigs/{read_sample}_to_{contig_sample}.{haplotype}.racon.rsrc'
    params:
        param_preset = lambda wildcards: config['minimap2_presets'][wildcards.read_sample]
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


rule quast_polished_assembly_statistics:
    input:
        contigs = 'output/assembly_polishing/racon_polished_contigs/{hapassm}.racon.ctg.fa',
        reference = 'references/' + REFNAME,
        genes = 'references/gencode.v30.basic.annotation.gff3.gz'
    output:
        pdf_report = 'output/assembly_analysis/quastlg/{hapassm}.racon/report.pdf',
        html_icarus = 'output/assembly_analysis/quastlg/{hapassm}.racon/icarus.html'
    threads: 16
    log: 'log/assembly_analysis/quastlg/{hapassm}.log'
    benchmark: 'run/assembly_analysis/quastlg/{hapassm}.rsrc'
    priority: 100
    params:
        output_dir = 'output/assembly_analysis/quastlg/{hapassm}.racon'
    run:
        exec = 'quast-lg.py --threads {threads}'
        exec += ' -r {input.reference}'
        exec += ' --features gene:{input.genes}'
        exec += ' --output-dir {params.output_dir}'
        exec += ' {input.contigs}'
        exec += ' &> {log}'
        shell(exec)

