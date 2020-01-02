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

Diagnostic plot of the BED output file can be produced as follows:

```R
library(SaaRclust)
library(BSgenome.Hsapiens.UCSC.hg38)

bedfile <- "bedfile with aligned contigs to the reference"
plt1 <- plotClusteredContigs(bedfile = bedfile, min.mapq = 10, bsgenome = BSgenome.Hsapiens.UCSC.hg38, report = 'clustering')
plt2 <- plotClusteredContigs(bedfile = bedfile, min.mapq = 10, bsgenome = BSgenome.Hsapiens.UCSC.hg38, report = 'ordering')
plt3 <- plotClusteredContigs(bedfile = bedfile, min.mapq = 10, bsgenome = BSgenome.Hsapiens.UCSC.hg38, report = 'orienting')

#To save the plots:
plot destination = "location and the file name where the plot should be saved"
ggsave(filename = <plot destination/clustering.pdf>, plot = plt1, width = 12, height = 6)
ggsave(filename = <plot destination/ordering.pdf>, plot = plt2, width = 12, height = 6)
ggsave(filename = <plot destination/orienting.pdf>, plot = plt3, width = 12, height = 6)
```