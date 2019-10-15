# Align Strand-seq reads to reference or custom assembly

Sent by David via e-mail on 2019-09-17

Relevance:
  - default alignment and preprocessing commands for strand-seq data
  - merge of mono- and dinucleotide fractions is data-specific

```
from collections import defaultdict

SAMPLE_DIR, SAMPLE, PLATE, LIBUM, NUCL = glob_wildcards(
    "rawData/{sample_dir}/{sample}_{plate}_{libum}_{nucl}_1.fastq.gz"
)

## Take unique ID for each sample
SAMPLES = sorted(set(SAMPLE))

rule all:
    input:
        bam=lambda wildcards: [
            "alignments/{}/{}.{}.{}.monodi.srt.mdup.bam.bai".format(
                sample_dir, sample, plate, libum
            ) for sample_dir, sample, plate, libum in zip(SAMPLE_DIR, SAMPLE, PLATE, LIBUM)
        ]

rule align_bwa:
    input:
        read1="rawData/{sample_dir}/{sample}_{plate}_{libum}_{nucl}_1.fastq.gz",
        read2="rawData/{sample_dir}/{sample}_{plate}_{libum}_{nucl}_2.fastq.gz",
        ref = config["reference"]
    output:
        temp("alignments/{sample_dir}/{nucl}/{sample}.{plate}.{libum}.{nucl}.bam")
    log:
        "log/{sample_dir}/{nucl}/{sample}.{plate}.{libum}.{nucl}.bam"
    threads:
        8
    params:
        rg="@RG\\tID:{sample}_{nucl}\\tPL:Illumina\\tSM:{sample}_{nucl}"
    shell:
        """
        bwa mem -t {threads} \
                -R '{params.rg}' \
                {input.ref} {input.read1} {input.read2} | samtools view -Sb -F 2304 - > {output} 2> {log}
        """

rule merge_mono_di:
    input:
        bam1="alignments/{sample_dir}/mono/{sample}.{plate}.{libum}.mono.bam",
        bam2="alignments/{sample_dir}/di/{sample}.{plate}.{libum}.di.bam"
    output:
        temp("alignments/{sample_dir}/{sample}.{plate}.{libum}.monodi.bam")
    log:
        "log/{sample_dir}/{sample}.{plate}.{libum}.monodi.bam.log"
    threads:
        8
    shell:
        "samtools merge -@ {threads} -O BAM {output} {input} 2> {log}"

rule sort_mono_di:
    input:
        "alignments/{sample_dir}/{sample}.{plate}.{libum}.monodi.bam"
    output:
        temp("alignments/{sample_dir}/{sample}.{plate}.{libum}.monodi.srt.bam")
    log:
        "log/{sample_dir}/{sample}.{plate}.{libum}.monodi.srt.bam.log"
    threads:
        8
    shell:
        "samtools sort -@ {threads} -O BAM {input} -o {output} 2> {log}"


rule markDuplicates:
    input:
        "alignments/{sample_dir}/{sample}.{plate}.{libum}.monodi.srt.bam"
    output:
        "alignments/{sample_dir}/{sample}.{plate}.{libum}.monodi.srt.mdup.bam"
    log:
        "alignments/{sample_dir}/{sample}.{plate}.{libum}.monodi.srt.mdup.bam.log"
    threads:
        8
    shell:
        "sambamba markdup -t {threads} {input} {output} 2> {log}"


rule index_bam:
    input:
        "alignments/{sample_dir}/{sample}.{plate}.{libum}.monodi.srt.mdup.bam"
    output:
        "alignments/{sample_dir}/{sample}.{plate}.{libum}.monodi.srt.mdup.bam.bai"
    threads:
        8
    shell:
        "samtools index -@ {threads} {input} {output}"

```
