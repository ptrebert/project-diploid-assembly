# Run and configure R tools

Sent by David on 2019-09-16

Relevance:
  - define several default parameters for R tools
  - reformatted for improved readability

```

rule run_saarclust:

## To run clustering
library(SaaRclust)
## BAMs with aligned Strand-seq reads against de novo assembly
bamfolder <- "/[...]/Saarclust_project/alignments_canu_HIFI/HG00733/"
## Location of the original assembly fasta file
assembly.fasta <- "/[...]/Saarclust_project/alignments_canu_HIFI/Assembly/<denovoAssembly.fasta>"
## Command to run SaaRclust scaffolding
scaffoldDenovoAssembly(
    bamfolder = bamfolder,
    outputfolder = "<outputfolder.path>/SaaRclust_results",
    store.data.obj = TRUE,
    reuse.data.obj = TRUE,
    assembly.fasta = assembly.fasta,
    concat.fasta = FALSE
)


rule run_breakpointr:

## Run breakpointR on clustered assembly to localize WC regions for phasing
inputfolder <- "BAM files with Strand-seq read aligned to the clustred assembly"
breakpointr(
    inputfolder = inputfolder,
    outputfolder = file.path(inputfolder, "BreakpointR_results"),
    windowsize = 500000,  # 500k
    binMethod = 'size',
    pairedEndReads = TRUE,
    pair2frgm = FALSE,
    min.mapq = 10,
    filtAlt = TRUE,
    background = 0.1,
    minReads = 50
)

## Export WC regions
exportRegions(
    datapath = file.path(inputfolder, "BreakpointR_results/data"),
    file = "HG00733_WCregion_clusteredAssembly.txt",
    collapseInversions = TRUE,
    collapseRegionSize = 5000000,  # 5M
    minRegionSize = 5000000,  # 5M
    state = "wc"
)


## To run StrandPhaseR manually
inputfolder <- "BAM files with Strand-seq read aligned to the clustred assembly"
strandPhaseR(
    inputfolder = inputfolder,
    positions = "path to SNV calls for a clustered assembly",
    WCregions = "file with exported wc regions by 'exportRegions' function",
    pairedEndReads = TRUE,
    min.mapq = 10,
    min.baseq = 20,
    num.iterations = 2,
    translateBases = TRUE,
    splitPhasedReads = TRUE,
    exportVCF ='HG00733'
)

rule run_strandphaser_per_chrom:
    input:
        wcregions    = "strand_states/{sample}/{windows}.{bpdens}/strandphaser_input.txt",
        snppositions = locate_snv_vcf,
        configfile   = "strand_states/{sample}/{windows}.{bpdens}/StrandPhaseR.{chrom}.config",
        strandphaser = "utils/R-packages/StrandPhaseR/R/StrandPhaseR",
        bamfolder    = "bam/{sample}/selected"
    output:
        "strand_states/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/StrandPhaseR_analysis.{chrom}/Phased/phased_haps.txt",
        "strand_states/{sample}/{windows}.{bpdens,selected_j[0-9\\.]+_s[0-9\\.]+_scedist[0-9\\.]+}/StrandPhaseR_analysis.{chrom}/VCFfiles/{chrom}_phased.vcf"
    log:
        "log/run_strandphaser_per_chrom/{sample}/{windows}.{bpdens}/{chrom}.log"
    shell:
        """
        Rscript utils/StrandPhaseR_pipeline.R \
                {input.bamfolder} \
                strand_states/{wildcards.sample}/{wildcards.windows}.{wildcards.bpdens}/StrandPhaseR_analysis.{wildcards.chrom} \
                {input.configfile} \
                {input.wcregions} \
                {input.snppositions} \
                $(pwd)/utils/R-packages/ \
                > {log} 2>&1
        """
```
