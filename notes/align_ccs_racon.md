# Align CCS reads to reference for Racon polishing

Sent by Aaron Wenger via e-mail on 2019-06-10

Relevance:
  - Racon accepts minimap's PAF output, but does not seem to perform polishing
  - Generally, for CCS reads, preset asm20 is suggested; not used here, see the following:

```
Use "-x map-pb" to minimap2 even when aligning CCS reads.
Racon will clip out segments that have no coverage, and map-pb does
a better job than `-x asm5` of avoiding alignment clipping at
quality dropouts in the draft assembly.
```

```
The parameters we use are:
$ minimap2 -a -x map-pb --eqx -m 5000 --secondary=no draft-asm.fa reads.fastq | samtools sort | samtools view -q 10 -F0x704 - > draft-asm.reads.sam
$ racon reads.fastq draft-asm.reads.sam draft-asm.fa -u > polished-asm.fa
```

