
import collections as col

REFKEY = config['use_ref_genome']
REFNAME = REFKEY + '.fa'
REFCHROM = config['chromosomes'][REFKEY]

SAMPLES = config['samples']

# TODO this needs to be reworked
REF_DOWNLOADS = ['GRCh38_decoy_hla.fa',
                 'GRCh38_giab_HG002_hconf-regions.bed',
                 'GRCh38_giab_HG002_hconf-variants.vcf.gz',
                 'GRCh38_giab_HG002_hconf-variants.vcf.gz.tbi']

# this needs to be sorted out when finalizing the pipeline
ruleorder: merge_phased_variants > prepare_trio_phased_variants
ruleorder: merge_chunk_part_read_files > merge_partitioned_read_files

include: 'smk_include/aux_utilities.smk'
include: 'smk_include/dl_read_data_chunks.smk'
include: 'smk_include/dl_read_data_parts.smk'

rule master:
    input:
        expand('references/{ref_file}',
                ref_file=REF_DOWNLOADS),

        # input validation
        expand('output/fastq_validation/{sample}.stats.pck',
                sample=SAMPLES),

        # produce sorted alignment files (against reference genome)
        expand('output/alignments/{sample}.sort.cram',
                sample=SAMPLES),

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
#         expand('output/consensus_haplotypes/{sample}.{project}.ont-ul.tag-{haplotype}{untagged}.ctg.fa',
#                sample=['HG002'], project=['pangen'],
#                haplotype=['h1', 'h2'], untagged=['', '-un']),
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


for record in config['ref_data']:
    if 'load' not in record:
        continue
    remote_url = record['load']['url']
    local_name = record['load']['name']
    local_path = 'references/' + local_name
    rule:
        output:
            local_path
        log: 'log/references/DL_{}.log'.format(local_name.rsplit('.', 1)[0])
        threads: 4
        message: 'Downloading reference file {}'.format(local_name)
        run:
            exec = 'aria2c -s {threads} -x {threads}'
            exec += ' -o {output}'
            exec += ' {}'.format(remote_url)
            exec += ' &> {log}'
            shell(exec)


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
    params:
        cache_path = os.path.join(os.getcwd(), config['ref_cache_folder']),
        ont_preset = lambda wildcards: 'map-ont' if 'ont' in wildcards.sample else '',
        pb_preset = lambda wildcards: 'map-pb' if 'pb' in wildcards.sample else ''
    threads: 48
    run:
        exec = 'REF_CACHE={params.cache_path}'
        exec += ' minimap2 -t {threads}'
        exec += ' -R \'@RG\\tID:1\\tSM:{wildcards.sample}\''
        exec += ' -a'  # output in SAM format
        exec += ' -x {params.ont_preset}{params.pb_preset}'
        exec += ' {input.reference}'
        exec += ' {input.read_data} 2> {log.minimap}'
        exec += ' | samtools view'
        exec += ' -T {input.reference}'
        exec += ' -C -'
        exec += ' > {output}'
        exec += ' 2> {log.samtools}'
        shell(exec)


# sorting takes ~1 hour / 50 GiB with 20 cores and 5GiB / core
# ATTENTION: current settings result in 20 * 20 GiB RAM => 400 GiB required!
rule sort_cram_alignments:
    input:
        'output/temp/alignments/{sample}.cram'
    output:
        'output/alignments/{sample}.sort.cram'
    log: 'log/aln_sort/{sample}.sort.log'
    threads: 20
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
        cram = 'output/sorted_aln/{sample}.sort.cram',
        crai = 'output/sorted_aln/{sample}.sort.cram.crai',
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


rule haplotag_reads:
    input:
        vcf = 'output/merged_phased_variants/{sample}.{project}.vcf.gz',
        tbi = 'output/merged_phased_variants/{sample}.{project}.vcf.gz.tbi',
        cram = 'output/sorted_aln/{sample}.{project}.ont-ul.sorted.cram',
        crai = 'output/sorted_aln/{sample}.{project}.ont-ul.sorted.cram.crai',
        ref = 'references/' + REFNAME,
    output:
        cram = 'output/tagged_aln/{sample}.{project}.ont-ul.sorted.tagged.cram',
        taglist = 'output/tagged_aln/{sample}.{project}.ont-ul.haplotag.tsv.gz'
    log: 'log/tagged_aln/{sample}.{project}.haplotag.log'
    conda:
        "environment/conda/wh_split.yml"
    shell:
        "whatshap --debug haplotag --output {output.cram} --reference {input.ref} --output-haplotag-list {output.taglist} {input.vcf} {input.cram} &> {log}"


rule haplosplit_fastq:
    input:
        fastq = 'input/read_data/aln_ready/{sample}.{project}.ont-ul.cmp.fastq.gz',
        taglist = 'output/tagged_aln/{sample}.{project}.ont-ul.haplotag.tsv.gz'
    output:
        haplo1 = 'output/haplosplit_reads/{sample}.{project}.ont-ul.tag-h1.fastq.gz',
        haplo2 = 'output/haplosplit_reads/{sample}.{project}.ont-ul.tag-h2.fastq.gz',
        untag = 'output/haplosplit_reads/{sample}.{project}.ont-ul.tag-un.fastq.gz'
    log:
        'log/haplosplit/{sample}.{project}.haplosplit.log'
    conda:
        "environment/conda/wh_split.yml"
    shell:
        "whatshap --debug split --pigz --output-h1 {output.haplo1} --output-h2 {output.haplo2} --output-untagged {output.untag} {input.fastq} {input.taglist} &> {log}"


rule assemble_haplotypes:
    input:
        fastq = 'output/haplosplit_reads/{sample}.{project}.ont-ul.tag-{haplotype}.fastq.gz',
    output:
        layout = 'output/assembled_haplotypes/{sample}.{project}.ont-ul.tag-{haplotype}.ctg.lay.gz',
    log: 'log/haplotype_assembly/{sample}.{project}.tag-{haplotype}.log'
    threads: 48
    run:
        exec = 'wtdbg2 -x ont'  # parameter preset for ONT
        exec += ' -i {input.fastq}'
        exec += ' -g3g -t {threads}'  # approx genome size
        exec += ' -o output/assembled_haplotypes/{wildcards.sample}.{wildcards.project}.ont-ul.tag-{wildcards.haplotype}'
        exec += ' &> {log}'
        shell(exec)


rule assemble_haplotypes_untagged:
    input:
        hap_fastq = 'output/haplosplit_reads/{sample}.{project}.ont-ul.tag-{haplotype}.fastq.gz',
        un_fastq = 'output/haplosplit_reads/{sample}.{project}.ont-ul.tag-un.fastq.gz',
    output:
        layout = 'output/assembled_haplotypes/{sample}.{project}.ont-ul.tag-{haplotype}-un.ctg.lay.gz',
    log: 'log/haplotype_assembly/{sample}.{project}.tag-{haplotype}-un.log'
    threads: 48
    run:
        exec = 'wtdbg2 -x ont'  # parameter preset for ONT
        exec += ' -i {input.hap_fastq} -i {input.un_fastq}'
        exec += ' -g3g -t {threads}'  # approx genome size
        exec += ' -o output/assembled_haplotypes/{wildcards.sample}.{wildcards.project}.ont-ul.tag-{wildcards.haplotype}-un'
        exec += ' &> {log}'
        shell(exec)


rule derive_assembly_consensus:
    input:
        ctg_layout = 'output/assembled_haplotypes/{hapassm}.ctg.lay.gz'
    output:
        'output/consensus_haplotypes/{hapassm}.ctg.fa'
    log: 'log/consensus_haplotypes/{hapassm}.log'
    threads: 48
    run:
        exec = 'wtpoa-cns -t {threads}'
        exec += ' -i {input.ctg_layout}'
        exec += ' -o {output} &> {log}'
        shell(exec)


rule compute_assembly_delta:
    input:
        contigs = 'output/consensus_haplotypes/{hapassm}.ctg.fa',
        reference = 'references/' + REFNAME
    output:
        delta = 'output/assembly_analysis/mummer/{hapassm}.delta'
    log: 'log/assembly_analysis/mummer/{hapassm}.log'
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
