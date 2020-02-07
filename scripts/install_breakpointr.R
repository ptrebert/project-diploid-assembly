#!/usr/bin/env Rscript

if (is.element('breakpointR', installed.packages()[,1])) {
    print('Removing previously installed version of breakpointR')
    remove.packages('breakpointR')
}

args = commandArgs(trailingOnly=TRUE)

git.commit = args[1]

if (is.na(as.numeric(git.commit))) {
    # means proper git tag

    # if dev version is installed, breakpointRdata
    # is not automatically included as a dependency,
    # so trigger setup manually
    devtools::install_git(
        "git://github.com/daewoooo/breakpointRdata.git",
        dependencies=FALSE,
        upgrade=FALSE
    )

    devtools::install_git(
        "git://github.com/daewoooo/breakpointR.git",
        ref = git.commit,
        dependencies=FALSE,
        upgrade=FALSE
    )
} else {
    BiocManager::install(
        c("breakpointR"),
        update=FALSE
    )
}

quit(save="no")