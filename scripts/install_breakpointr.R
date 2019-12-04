#!/usr/bin/env Rscript

if (is.element('breakpointR', installed.packages()[,1])) {
    print('Removing previously installed version of breakpointR')
    remove.packages('breakpointR')
}

BiocManager::install(
    c("breakpointR"),
    update=FALSE
)

quit(save="no")