#!/usr/bin/env Rscript

if (is.element('breakpointR', installed.packages()[,1])) {
    remove.packages('breakpointR')
}

BiocManager::install(c(
#    "Rsamtools",
    "breakpointR",
    update=FALSE
))

quit(save="no")