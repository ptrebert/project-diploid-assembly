#!/usr/bin/env Rscript

host.name = Sys.info()['nodename']

if (any(grepl('login', host.name, fixed=TRUE), grepl('submit', host.name, fixed=TRUE))) {
    stop('Running on cluster login/submit node - exiting')
}

suppressMessages(library(StrandPhaseR))

args = commandArgs(trailingOnly=TRUE)

commit.hash.expect = args[7]
if (!is.na(commit.hash.expect)) {
    pkg.desc = library(help=StrandPhaseR)
    pkg.desc.version = grep("RemoteRef", pkg.desc$info[[1]], ignore.case=TRUE, value=TRUE)
    pkg.desc.version = strsplit(pkg.desc.version, split=" +", fixed=FALSE)[[1]][2]
    install.match = grepl(commit.hash.expect, pkg.desc.version)
    if (!install.match) {
        stop(paste('Library StrandPhaseR version mismatch: installed >>>', pkg.desc.version, '<<< versus expected', commit.hash.expect, sep=' '))
    }
}

bam.folder = args[1]
breakpointr.data = args[2]
strandphaser.data = args[3]
strandphaser.vcf = args[4]
variant.calls = args[4]
clustered.assm = args[5]
output.folder = args[6]

correctInvertedRegionPhasing(
    input.bams = bam.folder,
    outputfolder = output.folder,
    inv.bed = NULL,
    recall.phased = TRUE,
    het.genotype = 'strict',
    snv.positions = variant.calls,
    breakpointR.data = breakpointr.data,
    strandphaseR.data = strandphaser.data,
    pairedEndReads = TRUE,
    min.mapq = 10,
    background = 0.1,
    vcfs.files = strandphaser.vcf,
    lookup.bp = 1000000
    ref.fasta = clustered.assm,
    assume.biallelic = TRUE
)

warnings()

quit(save='no')