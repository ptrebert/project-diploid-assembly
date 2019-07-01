

rule build_flye_assembly_master:
    input:
        expand('output/haplotype_assembly/consensus/{sample}.{haplotype}.flye/{sample}.{haplotype}.ctg.fa',
                sample=['HG00514.sequel2.pb-clr-25'], haplotype=['tag-h1-un', 'tag-h2-un'])


onsuccess:
    body_text = "Nothing to report\n{log}"
    if config['notify']:
        shell('mail -s "[Snakemake] flye assembly - success" {} <<< "{}"'.format(config['notify_email'], body_text))

onerror:
    if config['notify']:
        shell('mail -s "[Snakemake] flye assembly - ERRROR" {} < {{log}}'.format(config['notify_email']))


rule build_flye_haplotype_assembly:
    input:
        'output/haplotype_partitioning/splitting/{sample}.{haplotype}.fastq.gz'
    output:
        'output/haplotype_assembly/consensus/{sample}.{haplotype}.flye/{sample}.{haplotype}.ctg.fa'
    log: 'log/haplotype_assembly/consensus/{sample}.{haplotype}.flye.log'
    benchmark: 'run/haplotype_assembly/consensus/{sample}.{haplotype}.flye.rsrc'
    threads: 48
    params:
        readset = lambda wildcards: config['flye_presets'][wildcards.sample],
        output_dir = 'output/haplotype_assembly/consensus/{sample}.{haplotype}.flye'
    conda:
        "../environment/conda/conda_pbtools.yml"
    shell:
        'flye --{params.readset} {input} --genome-size 3g --threads {threads} --iterations 2 --debug --out-dir {params.output_dir} &> {log}'