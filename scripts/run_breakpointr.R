#!/usr/bin/env Rscript

suppressMessages(library(breakpointR))

args = commandArgs(trailingOnly=TRUE)

bam.folder = args[1]
config.file = args[2]
output.folder = args[3]
num.cpu = args[4]
output.file = args[5]

breakpointr(
    inputfolder=bam.folder,
    outputfolder=output.folder,
    configfile=config.file,
    numCPU=num.cpu
)

exportRegions(
    datapath=file.path(output.folder, "data"),
    file=output.file,
    collapseInversions=TRUE,
    collapseRegionSize=5000000,
    minRegionSize=5000000,
    state="wc"
)

quit(save='no')