# Snakefile: minimap contig to reference assembly alignment

Sent via e-mail by David on 2019-10-08

Relevance:
  - minimap2 parameters to get reasonable alignments between de novo and reference
  - reference (`ref` below) refers here to, e.g., GRCh38


```
    ## Snakefile to align denovo contigs to the reference genome

    ## Set config file
    configfile: "Snake.config.json"

    FASTA, = glob_wildcards("clustered_assembly/{fasta}.fasta")

    rule all:
        input:
            #expand("alignments/{fasta}.bed", fasta=FASTA)
            "alignments/HG00733_sra_pbsq1-clr_sqa_clustered_v2.bed"

    rule align_fasta:
        input:
            fasta = "clustered_assembly/{fasta}.fasta",
            ref = config["reference"]
        output:
            "alignments/{fasta}.bam"
        log:
            "log/{fasta}.bam.log"
        threads:
            8
        shell:
            #"minimap2 -ax asm20 --eqx -r 20000 -s 30000 -t {threads} {input.ref} {input.fasta} | samtools view -F 260 -b - | samtools sort - > {output} 2> {log}"
            "minimap2 --secondary=no --eqx -Y -ax asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 -t {threads} {input.ref} {input.fasta} | samtools view -F 260 -b - | samtools sort - > {output} 2> {log}"

    rule bam2bed:
        input:
            "alignments/{fasta}.bam"
        output:
            "alignments/{fasta}.bed"
        log:
            "log/{fasta}.bed.log"
        shell:
            "bedtools bamtobed -i {input} > {output} 2> {log}"

```
