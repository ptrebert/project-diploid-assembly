#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

git.commit = args[1]

BiocManager::install(c(
    "GenomicRanges",  # 1.34.0
    "IRanges",  # 2.16.0
    "Rsamtools",  # 1.34.1
    "GenomicAlignments",  # 1.18.1
    "GenomeInfoDb",  # 1.18.2
    "GenomeInfoDbData",  # 1.2.0,
    update=FALSE
))

devtools::install_git("git://github.com/daewoooo/SaaRclust.git", ref = git.commit)

quit(save="no")