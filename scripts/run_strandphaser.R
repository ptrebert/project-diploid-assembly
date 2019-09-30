#!/usr/bin/env Rscript

suppressMessages(library(StrandPhaseR))

args = commandArgs(trailingOnly=TRUE)

bam.folder = args[1]
config.file = args[2]
variant.calls = args[3]
wc.regions = args[4]
output.folder = args[5]
sample.individual = args[6]

strandPhaseR(
    inputfolder=bam.folder,
    configfile=config.file,
    outputfolder=output.folder,
    positions=variant.calls,
    WCregions=wc.regions,
    exportVCF=sample.individual
)

quit(save='no')