#!/usr/bin/env Rscript

if (is.element('breakpointR', installed.packages()[,1])) {
    remove.packages('breakpointR')
}

BiocManager::install(c(
    "Rsamtools",
    "breakpointR"
))

quit(save="no")