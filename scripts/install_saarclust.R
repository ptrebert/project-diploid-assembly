#!/usr/bin/env Rscript

BiocManager::install(c(
    "GenomicRanges",  # 1.34.0
    "IRanges",  # 2.16.0
    "Rsamtools",  # 1.34.1
    "GenomicAlignments",  # 1.18.1
    "GenomeInfoDb",  # 1.18.2
    "GenomeInfoDbData"  # 1.2.0
))

devtools::install_git("git://github.com/daewoooo/SaaRclust.git", ref = "devel")
