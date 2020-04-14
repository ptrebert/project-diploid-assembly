#!/usr/bin/env Rscript

suppressMessages(library(SaaRclust))

args = commandArgs(trailingOnly=TRUE)

scaffoldDenovoAssembly(
    configfile = args[1],
    bamfolder = args[2],
    outputfolder = args[3]
)

warnings()

quit(save='no')